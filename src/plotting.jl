using Printf
@userplot Specificationplot

"""
    specificationplot([S,CS,T], [WS,WU,WT])

This function visualizes the control synthesis using the hInf_synthesize with
the three weighting functions {WS(jω), WU(jω), WT(jω)} inverted and scaled by γ,
against the corresponding transfer fucntions {S(jω), C(jω)S(jω), T(jω)}, to
verify visually that the specifications are met. This may be run using both MIMO
and SISO systems.
"""
specificationplot
@recipe function specificationplot(
    p::Specificationplot;
    wint = (-3, 5),
    wnum = 101,
    hz = true,
    nsigma=typemax(Int),
    s_labels = [
        "σ(S(jω))",
        "σ(C(jω)S(jω))",
        "σ(T(jω))",
    ],
    w_labels = [
        "γ σ(Wₛ(jω)⁻¹)",
        "γ σ(Wᵤ(jω)⁻¹)",
        "γ σ(Wₜ(jω)⁻¹)",
    ],
    colors = [:red, :blue, :green],
)

    sensitivityfunctions = p.args[1]
    if length(p.args) >= 3
        weightfunctions, γ = p.args[2:3]
    end

    title --> "Specification sigma plot"
    xguide --> "Frequency ($(hz ? "Hz" : "rad/s"))", yguide --> "Singular Values $_PlotScaleStr"

    w = [10^i for i in range(wint[1], stop = wint[2], length = wnum)]

    # Plot the sensitivity functions
    for (index, G) in enumerate(sensitivityfunctions)
        if G isa Number || G isa LTISystem
            singval = sigma(ss(G), w)[1]
            if ControlSystems._PlotScale == "dB"
                singval = 20 * log10.(singval)
            end
            for i = 1:min(size(singval, 2), nsigma)
                @series begin
                    xscale --> :log10
                    yscale --> ControlSystems._PlotScaleFunc
                    linestyle --> :solid
                    color --> colors[mod(index - 1, 3)+1]
                    label --> (i == 1 ? s_labels[mod(index - 1, 3)+1] : "")
                    w ./ (hz ? 2pi : 1), singval[:, i]
                end
            end
        end
    end
    if length(p.args) >= 3

        ## Plot the weight functions
        for (index, W) in enumerate(weightfunctions)
            if W isa Number
                W = ss(float(W))
            elseif W == []
                continue
            end
            singval = sigma(γ / W, w)[1]
            if ControlSystems._PlotScale == "dB"
                singval = 20 * log10.(singval)
            end
            for i = 1:size(singval, 2)
                weightlabel = (i == 1 ? w_labels[mod(index - 1, 3)+1] : "")
                @series begin
                    xscale --> :log10
                    yscale --> ControlSystems._PlotScaleFunc
                    linestyle --> :dash
                    color --> colors[mod(index - 1, 3)+1]
                    linewidth --> 2
                    label --> (i == 1 ? w_labels[mod(index - 1, 3)+1] : "")
                    w ./ (hz ? 2pi : 1), singval[:, i]
                end
            end
        end
    end
end

# Case where a single sensitivity function (for instance the closed loop TF from
# disturbance to output) and the gain γ
specificationplot(sens::T, γ::Number; kwargs...) where {T<:LTISystem} =
    specificationplot([sens], [1], γ; kwargs...)




@userplot MvNyquistplot
"""
    fig = nyquistplot(sys, w;  unit_circle=false, hz = false, kwargs...)

Create a Nyquist plot of the `LTISystem`. A frequency vector `w` must be
provided.

- `unit_circle`: if the unit circle should be displayed
If `hz=true`, the hover information will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to plot.
"""
nyquistplot
@recipe function nyquistplot(p::MvNyquistplot;  unit_circle=false, hz=false)
    sys, w = p.args
    nw = length(w)
    framestyle --> :zerolines
    θ = range(0, stop=2π, length=100)
    S, C = sin.(θ), cos.(θ)
    L = freqresp(sys, w)
    dets = map(axes(L, 1)) do i
        det(I + L[i,:,:])
    end
    dets = vec(dets)
    redata = real.(dets)
    imdata = imag.(dets)
    if eltype(imdata) <: AbstractFloat
        ylims --> (clamp(minimum(imdata), -10, -1), clamp(maximum(imdata), 1, 10))
        xlims --> (clamp(minimum(redata), -10, -1), clamp(maximum(redata), 1, 10))
    end
    title --> "Mutivariable Nyquist plot"

    @series begin
        hover --> [hz ? Printf.@sprintf("f = %.3f", w/2π) : Printf.@sprintf("ω = %.3f", w) for w in w]
        (redata, imdata)
    end                
    @series begin # Mark the critical point
        primary := false
        markershape := :xcross
        seriescolor := :red
        markersize := 5
        seriesstyle := :scatter
        [0], [0]
    end             
    if unit_circle 
        @series begin
            primary := false
            linestyle := :dash
            linecolor := :gray
            seriestype := :path
            markershape := :none
            (C, S)
        end
    end
end


@userplot Muplot

"""
    muplot(sys, args...; hz=false)
    muplot(LTISystem[sys1, sys2...], args...; hz=false)

Plot the structured singular values of the frequency response of the `LTISystem`(s). This plot is similar to `sigmaplot`, but scales the loop-transfer function to minimize the maximum singular value. Only applicable to square systems.
A frequency vector `w` can be optionally provided.

If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to Plots.plot.
"""
muplot
@recipe function muplot(p::Muplot; hz=false)
    systems, w = ControlSystems._processfreqplot(Val{:sigma}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    nw = length(w)
    title --> "Structured singular value plot (μ)"
    xguide --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
    yguide --> "μ Singular Values $(ControlSystems._PlotScaleStr)"
    @views for (si, s) in enumerate(systems)
        M = freqresp(s, w)
        M = permutedims(M, (2,3,1))
        mu, D = structured_singular_value(M, scalings=true, tol=1e-6, dynamic=true)
        Di = Diagonal(D)
        sv = map(axes(M, 3)) do i
            svdvals(Di\M[:,:,i]*Di)
        end

        sv = reduce(hcat, sv)'
        if ControlSystems._PlotScale == "dB"
            @. sv = 20*log10(sv)
        end
        for i in 1:size(sv, 2)
            @series begin
                xscale --> :log10
                yscale --> ControlSystems._PlotScaleFunc
                color --> si
                ws, sv[:, i]
            end
        end
    end
end
