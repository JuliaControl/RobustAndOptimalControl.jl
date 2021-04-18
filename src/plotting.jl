
@userplot Specificationplot
"""`specificationplot(Any[S,CS,T], Any[WS,WU,WT], args...)`

This function visualizes the control synthesis using the hInf_synthesize with
the three weighting functions {WS(jω), WU(jω), WT(jω)} inverted and scaled by γ,
against the the corresponding transfer fucntions {S(jω), C(jω)S(jω), T(jω)}, to
verify visually that the specifications are met. This may be run using both MIMO
and SISO systems.

`kwargs` is sent as argument to Plots.plot."""
specificationplot
@recipe function specificationplot(
    p::Specificationplot;
    wint = (-3, 5),
    wnum = 101,
    s_labels = [
        "\$\\sigma (S(j\\omega))\$",
        "\$\\sigma (C(j\\omega)S(j\\omega))\$",
        "\$\\sigma (T(j\\omega))\$",
    ],
    w_labels = [
        "\$\\gamma \\sigma (W_S(j\\omega)^{-1})\$",
        "\$\\gamma \\sigma (W_U(j\\omega)^{-1})\$",
        "\$\\gamma \\sigma (W_T(j\\omega)^{-1})\$",
    ],
    colors = [:red, :blue, :green],
)

    sensitivityfunctions, weightfunctions, gamma = p.args[1:3]

    title --> "Specification sigma plot"
    xguide --> "Frequency (rad/s)", yguide --> "Singular Values $_PlotScaleStr"

    w = [10^i for i in range(wint[1], stop = wint[2], length = wnum)]

    # Plot the sensitivity functions
    for (index, G) in enumerate(sensitivityfunctions)
        if G isa Number || G isa LTISystem
            singval = sigma(ss(G), w)[1]
            if ControlSystems._PlotScale == "dB"
                singval = 20 * log10.(singval)
            end
            for i = 1:size(singval, 2)
                @series begin
                    xscale --> :log10
                    yscale --> ControlSystems._PlotScaleFunc
                    linestyle --> :dash
                    linecolor --> colors[mod(index - 1, 3)+1]
                    label --> (i == 1 ? s_labels[mod(index - 1, 3)+1] : "")
                    w, singval[:, i]
                end
            end
        end
    end

    ## Plot the weight functions
    for (index, W) in enumerate(weightfunctions)
        if W isa Number
            W = ss(W)
        end
        if W isa LTISystem
            if size(W) != (1, 1)
                error(
                    ErrorException(
                        "We can currently only handle SISO weight funcitions in the visualization",
                    ),
                )
            end
            singval = sigma(gamma / tf(W), w)[1]
            if ControlSystems._PlotScale == "dB"
                singval = 20 * log10.(singval)
            end
            for i = 1:size(singval, 2)
                weightlabel = (i == 1 ? w_labels[mod(index - 1, 3)+1] : "")
                @series begin
                    xscale --> :log10
                    yscale --> ControlSystems._PlotScaleFunc
                    linestyle --> :dot
                    linecolor --> colors[mod(index - 1, 3)+1]
                    linewidth --> 2
                    label --> (i == 1 ? w_labels[mod(index - 1, 3)+1] : "")
                    w, singval[:, i]
                end
            end
        end
    end
end

specificationplot(
    sens::Vector{T},
    weight::Vector{T},
    gamma::Number;
    kwargs...,
) where {T<:LTISystem} = specificationplot(sens, weight, gamma; kwargs...)

# Case where a ingle sensitivity function (for instance the closed loop TF from
# disturbance to output) and the gain γ
specificationplot(sens::T, gamma::Number; kwargs...) where {T<:LTISystem} =
    specificationplot([sens], [1], gamma; kwargs...)

