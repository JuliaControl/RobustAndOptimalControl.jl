"""
    Diskmargin

The notation follows "An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet

# Fields:
`α`: The disk margin
`ω0`: The worst-case frequency
`f0`: The destabilizing perturbation `f0` is a complex number with simultaneous gain and phase variation. This critical perturbation causes an instability with closed-loop pole on the imaginary axis at the critical frequency ω0 
`δ0`: The uncertain element generating f0.
`γmin`: The lower real-axis intercept of the disk (analogous to classical lower gain margin).
`γmax`: The upper real-axis intercept of the disk (analogous to classical upper gain margin).
`ϕm`: is analogous to the classical phase margin.
`σ`: The skew parameter that was used to calculate the margin

Note, `γmax` and `ϕm` are in smaller than the classical gain and phase margins sicne the classical margins do not consider simultaneous perturbations in gain and phase. 

The "disk" margin becomes a half plane for `α = 2` and an inverted circle for `α > 2`. In this case, the upper gain margin is infinite. See the paper for more details, in particular figure 6.
"""
struct Diskmargin
    α
    ω0
    f0
    δ0
    γmin
    γmax
    σ
    ϕm
    L
end

function Diskmargin(α, σ=0; ω0=mising, f0=missing, δ0=missing, L=nothing)
    d = Disk(; α, σ)
    γmin = d.γmin
    γmax = d.γmax
    ϕm = d.ϕm
    Diskmargin(α, ω0, f0, δ0, γmin, γmax, σ, ϕm, L)
end


function Base.getproperty(dm::Diskmargin, s::Symbol)
    s ∈ fieldnames(typeof(dm)) && return getfield(dm, s)
    if s ∈ (:margin, :alpha, :diskmargin)
        return dm.α
    elseif s ∈ (:gm, :gainmargin)
        return [dm.γmin, dm.γmax]
    elseif s === :phasemargin
        return dm.ϕm
    else
        throw(ArgumentError("$(typeof(dm)) has no property named $s"))
    end
end

Base.propertynames(dm::Diskmargin) = (fieldnames(typeof(dm))..., :margin, :gainmargin, :phasemargin)


function Base.show(io::IO, dm::Diskmargin)
    println(io, "Disk margin with:")
    println(io, "Margin: ", dm.α)
    println(io, "Frequency: ", dm.ω0, " rad/s,  ", dm.ω0/(2π), " Hz")
    if dm.γmax < dm.γmin # In this case, we have an "inverted circle"
        println(io, "Gain margins: [$(dm.γmin), Inf]")
    else
        println(io, "Gain margins: [$(dm.γmin), $(dm.γmax)]")
    end
    println(io, "Phase margin: ", dm.ϕm)
    delaymarg = π/180 * dm.ϕm / dm.ω0
    print(io, "Delay margin: ", delaymarg, " s")
    if isdiscrete(dm.L)
        samples = delaymarg / dm.L.Ts
        println(io, ",  ", isfinite(samples) ? floor(Int, samples) : samples, " samples")
    else
        println(io)
    end
    println(io, "Skew: ", dm.σ)
    println(io, "Worst-case perturbation: ", dm.f0)
end

"""
    Disk

Represents a perturbation disc in the complex plane. `Disk(0.5, 2)` represents all perturbations in the circle centered at 1.25 with radius 0.75, or in other words, a gain margin of 2 and a phase margin of 36.9 degrees.

A disk can be converted to a Nyquist exclusion disk by `nyquist(disk)` and plotted using `plot(disk)`.

# Arguments:
- `γmin`: Lower intercept
- `γmax`: Upper intercept
- `c`: Center
- `r`: Radius
- `ϕm`: Angle of tangent line through origin.

If γmax < γmin the disk is inverted.
See [`diskmargin`](@ref) for disk margin computations. 
"""
struct Disk{T}
    γmin::T
    γmax::T
    c::T
    r::T
    ϕm::T
end

Disk(args...) = Disk(promote(args...)...)

center_radius(γmin, γmax) = 1/2 * (γmax + γmin), 1/2 * (γmax - γmin)

function Disk(γmin, γmax)
    c, r = center_radius(γmin, γmax)
    Disk(γmin, γmax, c, r)
end

function Disk(γmin, γmax, c, r)
    if !isfinite(γmax)
        ϕm = 90.0
    else
        ϕm = (1 + γmin*γmax) / (γmin + γmax)
        ϕm = ϕm >= 1 ? 0 : ϕm <= -1 ? Inf : rad2deg(acos(ϕm))
    end
    Disk(γmin, γmax, c, r, ϕm)
end

function Disk(; α, σ)
    γmin = (2 - α*(1-σ)) / (2 + α*(1+σ))
    γmax = (2 + α*(1-σ)) / (2 - α*(1+σ))
    Disk(γmin, γmax)
end

Disk(dm::Diskmargin) = Disk(dm.γmin, dm.γmax)

"""
    nyquist(d::Disk)

Transform a `Disk` representing a diskmargin to a exclusion disk in the Nyquist plane. This can be useful for visualizing a diskmargin in the Nyquist plane.
"""
ControlSystemsBase.nyquist(d::Disk) = Disk(-inv(d.γmin), -inv(d.γmax)) # translate the disk to a nyquist exclusion disk

"""
    diskmargin(L, σ = 0)
    diskmargin(L, σ::Real, ω)

Calculate the disk margin of LTI system `L`. `L` is supposed to be a loop-transfer function, i.e., it should be square. If `L = PC` the disk margin for output perturbations is computed, whereas if `L = CP`, input perturbations are considered. If the method `diskmargin(P, C, args...)` is used, both are computed. Note, if `L` is MIMO, a simultaneous margin is computed, see [`loop_diskmargin`](@ref) for single loop margins of MIMO systems.

The implementation and notation follows
["An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet](https://arxiv.org/abs/2003.04771).


The margins are aviable as fields of the returned objects, see [`Diskmargin`](@ref).

# Arguments:
- `L`: A loop-transfer function.
- `σ`: If little is known about the distribution of gain variations then σ = 0 is a reasonable choice as it allows for a gain increase or decrease by the same relative amount. *The choice σ < 0* is justified if the gain can decrease by a larger factor than it can increase. Similarly, *the choice σ > 0* is justified when the gain can increase by a larger factor than it can decrease. *If σ = −1* then the disk margin condition is αmax = inv(MT). This margin is related to the robust stability condition for models with multiplicative uncertainty of the form P (1 + δ). If σ = +1 then the disk margin condition is αmax = inv(MS)
- `kwargs`: Are sent to the [`hinfnorm`](@ref) calculation
- `ω`: If a vector of frequencies is supplied, the frequency-dependent disk margin will be computed, see example below.

# Example: 
```
L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
plot(dm) # Plot the disk margin to illustrate maximum allowed simultaneous gain and phase variations.

nyquistplot(L)
plot!(dm, nyquist=true) # plot a nyquist exclusion disk. The Nyquist curve will be tangent to this disk at `dm.ω0`
nyquistplot!(dm.f0*L) # If we perturb the system with the worst-case perturbation `f0`, the curve will pass through the critical point -1.

## Frequency-dependent margin
w = exp10.(LinRange(-2, 2, 500))
dms = diskmargin(L, 0, w)
plot(dms; lower=true, phase=true)
```

# Example: relation to Ms and Mt
```
Ms, wMs = hinfnorm(input_sensitivity(P, C)) # Input Ms
dm = diskmargin(C*P, 1) # Input diskmargin, skew = +1
isapprox(Ms/(Ms-1), dm.gainmargin[2], rtol=1e-2) # Guaranteed gain margin based on Ms
isapprox(inv(Ms), dm.margin, rtol=1e-2)
isapprox(dm.ω0, wMs, rtol=1e-1)


Mt, wMt = hinfnorm(input_comp_sensitivity(P, C)) # Input Mt
dm = diskmargin(C*P, -1) # Input diskmargin, skew = -1
isapprox(inv(Mt), dm.margin, rtol=1e-2)
isapprox(dm.ω0, wMt, rtol=1e-1)
```

See also [`ncfmargin`](@ref) and [`loop_diskmargin`](@ref).
"""
function diskmargin(L::LTISystem, σ::Real=0; l=1e-3, u=1e3, kwargs...)
    L isa DelayLtiSystem && @warn "To compute the diskmargin of delay systems, consider approximating the delay with a Pade-approximation (`pade`) or discretize the system (`c2d`)."
    issiso(L) || return sim_diskmargin(L, σ, l, u)
    M = feedback(1, L) + (σ-1)/2
    n,ω0 = hinfnorm2(M; kwargs...)
    diskmargin(L, σ, ω0)
end

"""
    diskmargin(L::LTISystem, σ::Real, ω)

Calculate the diskmargin at a particular frequency or vector of frequencies. If `ω` is a vector, you get a frequency-dependent diskmargin plot if you plot the returned value.
See also [`ncfmargin`](@ref).
"""
diskmargin(L::LTISystem, σ::Real, ω::AbstractArray) = map(w->diskmargin(L, σ, w), ω)

function diskmargin(L::LTISystem, σ::Real, ω0::Real)
    L isa DelayLtiSystem && @warn "To compute the diskmargin of delay systems, consider approximating the delay with a Pade-approximation (`pade`) or discretize the system (`c2d`)."
    issiso(L) || return sim_diskmargin(L, σ, [ω0])[]
    M = feedback(1, L) + (σ-1)/2
    freq = isdiscrete(L) ? cis(ω0*L.Ts) : complex(0, ω0)
    Sω = evalfr(M, freq)[]
    αmax = 1/abs(Sω)
    δ0 = inv(Sω)
    dp = Disk(; α = αmax, σ)
    if δ0 == 2/(σ+1)  # S = 1, L = 0
        Diskmargin(αmax, ω0, Inf, δ0, 0, Inf, σ, dp.ϕm, L)
    else
        f0 = (2 + δ0*(1-σ)) / (2 - δ0*(1+σ))
        Diskmargin(αmax, ω0, f0, δ0, dp.γmin, dp.γmax, σ, dp.ϕm, L)
    end
end



## Plotting ====================================================================


γϕcurve(dm::Diskmargin; kwargs...) = γϕcurve(dm.α, dm.σ; kwargs...)

"""
    γϕcurve(α, σ; N = 200)

Internal function. Get the curve of extremal stable gain and phase perturbations. This function is called with a single diskmargin object is plotted.

- `N`: Number of points on the curve.
"""
function γϕcurve(α, σ; N = 200)
    θ = LinRange(0, π, N)
    f = @. (2 - α*cis(θ)*(1-σ)) / (2 + α*cis(θ)*(1+σ))
    @. (abs(f), abs(rad2deg(angle(f))))
end

@recipe function plotdisk(d::Disk; nyquist = false)
    θ = LinRange(0, 2pi, 200)
    re, im = @. cos(θ), sin(θ)
    if nyquist
        d = ControlSystemsBase.nyquist(d)
    end
    c,r = d.c, d.r
    @series begin 
        fill --> true
        fillalpha --> 0.5
        @. r*re+c, r*im
    end
end

@recipe function plotdiskmargin(dm::Diskmargin; nyquist=false)
    if nyquist
        @series begin
            label --> "σ = $(dm.σ)"
            ControlSystemsBase.nyquist(Disk(dm))
        end
    else
        γ, ϕ = γϕcurve(dm)
        @series begin
            title --> "Stable region for combined gain and phase variation"
            xguide --> "Gain variation"
            yguide --> "Phase variation"
            label --> "σ = $(dm.σ)"
            fill --> true
            fillalpha --> 0.5
            γ, ϕ
        end
    end
end

@recipe function plotdiskmargin(dm::AbstractVector{<:Diskmargin}; lower=true, phase=true, gain=true)
    gain || phase || error("Must plot either gain or phase")
    w = [dm.ω0 for dm in dm]
    layout --> (phase + gain, 1)
    link --> :x
    gma = getfield.(dm, :γmax)
    gmi = getfield.(dm, :γmin)
    # ylims --> (0, min(10, maximum(gma)))
    data = (lower ? [gmi gma] : gma)
    data = max.(0, data)
    replace!(data, 0 => -Inf)
    if gain
        @series begin
            subplot --> 1
            label --> (lower ? ["Lower" "Upper"] : "Upper")
            # xguide --> "Frequency"
            xscale --> :log10
            yscale --> :log10
            w, data
        end
        replace!(data, -Inf => Inf)
        m,i = findmin(data[:, end])
        @series begin
            # subplot --> 1
            title --> "Gain margin"
            primary := true
            seriestype := :scatter
            seriescolor := :red
            label := string(round(m, sigdigits=3))
            [w[i]], [m]
        end
        @series begin
            primary := false
            seriestype := :hline
            linecolor := :black
            [1]
        end
    end
    if phase
        ϕm = getfield.(dm, :ϕm)
        m,i = findmin(ϕm)
        @series begin
            subplot --> (phase + gain)
            primary := true
            seriestype := :scatter
            seriescolor := :red
            label := string(round(m, sigdigits=3))
            [w[i]], [m]
        end
        @series begin
            subplot --> (phase + gain)
            title --> "Phase margin (deg)"
            xguide --> "Frequency"
            xscale --> :log10
            label --> ""
            w, ϕm
        end
    end
end

@recipe function plotdiskmarginvector(dm::AbstractVector{<:AbstractVector{Diskmargin}})
    reduce(hcat, dm)
end


@recipe function plotdiskmarginmatrix(dm::AbstractMatrix{Diskmargin}; lower=true, phase=true)
    w = getfield.(dm[:,1], :ω0)
    ny = size(dm, 2)
    length(w) == size(dm, 1) || throw(ArgumentError("Frequency vector and diskmargin vector must have the same lengths."))
    layout --> (phase ? 2 : 1, 1)
    link --> :x
    neg2inf(x) = x <= 0 ? Inf : x
    for i = 1:ny
        is = ny == 1 ? "" : "$i"
        @series begin
            subplot --> 1
            title --> "Gain margin"
            # label --> permutedims([["Lower $i" for i = 1:ny]; ["Upper $i" for i = 1:ny]])
            label --> (lower ? ["Upper"*is "Lower"*is] : "Upper"*is)
            # xguide --> "Frequency"
            xscale --> :log10
            yscale --> :log10
            gma = getfield.(dm[:, i], :γmax) .|> neg2inf
            gmi = getfield.(dm[:, i], :γmin) .|> neg2inf
            replace!(gmi, 0 => -Inf)
            # ylims --> (0, min(10, maximum(gma)))
            w, (lower ? [gmi gma] : gma)
        end
        phase || continue
        @series begin
            subplot --> 2
            title --> "Phase margin"
            xguide --> "Frequency"
            xscale --> :log10
            label --> is#permutedims(["$i" for i = 1:ny])
            ϕm = getfield.(dm[:, i], :ϕm)
            w, ϕm
        end
    end
end

@userplot GainPhasePlot

@recipe function f(gp::GainPhasePlot)
    if length(gp.args) == 3
        P, re, im = gp.args
    elseif length(gp.args) == 1
        P = gp.args[1]
        re = LinRange(-1.5, 10, 300)
        im = LinRange(-3, 3, 300)
    else
        error("gainphaseplot takes either 1 (sys) or 3 arguments (sys, re, im).")
    end

    res = map(Iterators.product(re,im)) do (r,i)
        isstable(feedback(P*complex(r,i)))
    end
    colorbar --> false
    aspect_ratio --> 1
    @series begin
        color --> [:red, :green]
        seriesalpha --> 0.5
        seriestype := :heatmap
        re, im, res'
    end
    θ = LinRange(0, 2π, 100)
    S = sin.(θ)
    C = cos.(θ)
    framestyle := :zerolines
    @series begin
        primary := false
        linestyle := :dash
        seriescolor := :black
        C, S
    end
    r = [0, 1.5]
    c = [:red, :orange, :green]
    for (i, θ) = enumerate([30, 45, 60])
        points = r*cis(deg2rad(θ))
        @series begin
            linestyle := :dash
            linewidth --> 2
            xguide --> "Re"
            yguide --> "Im"
            seriescolor --> c[i]
            label --> string(θ)*"°"
            real(points), imag(points)
        end
    end
end

"""
    gainphaseplot(P)
    gainphaseplot(P, re, im)

Plot complex perturbantions to the plant `P` and indicate whether or not the closed-loop system is stable. The diskmargin is the largest disk that can be fit inside the green region that only contains stable variations.
"""
gainphaseplot


"""
    passivity_index(P; kwargs...)

Return
```math
γ = \\begin{Vmatrix}
(I-P)(I+P)^{-1}
\\end{Vmatrix}_∞
```
If ``γ ≤ 1``, the system is passive. If the system has unstable zeros, ``γ = ∞``

The negative feedback interconnection of two systems with passivity indices γ₁ and γ₂ is stable if ``γ₁γ₂ < 1``.

A passive system has a Nyquist curve that lies completely in the right half plane, and satisfies the following inequality (dissipation of energy)
```math
\\int_0^T y^T u dt > 0 ∀ T
```
The negative feedback-interconnection of two passive systems is stable and  parallel connections of two passive systems as well as the inverse of a passive system are also passive. A passive controller will thus always yeild a stable feedback loop for a passive system. A series connection of two passive systems *is not* always passive.

See also [`ispassive`](@ref), [`passivityplot`](@ref).
"""
function passivity_index(P; kwargs...)
    P.ny == P.nu || throw(ArgumentError("passivity_index only defined for square systems"))
    G = (ss(I(P.ny), P.timeevol)-P)*feedback(1, P)
    hinfnorm2(G; kwargs...)[1]
end

"""
    ispassive(P; kwargs...)

Determine if square system `P` is passive, i.e., ``P(s) + Pᴴ(s) > 0``.

A passive system has a Nyquist curve that lies completely in the right half plane, and satisfies the following inequality (dissipation of energy)
```math
\\int_0^T y^T u dt > 0 ∀ T
```
The negative feedback-interconnection of two passive systems is stable and  parallel connections of two passive systems as well as the inverse of a passive system are also passive. A passive controller will thus always yeild a stable feedback loop for a passive system. A series connection of two passive systems *is not* always passive.

See also [`passivityplot`](@ref), [`passivity_index`](@ref).
"""
ispassive(G; kwargs...) = passivity_index(G; kwargs...) <= 1

@userplot Passivityplot
"""
    passivityplot(sys, args...; hz=false)
    passivityplot(LTISystem[sys1, sys2...], args...; hz=false)

Plot the passivity index of a `LTISystem`(s). The system is passive for frequencies where the index is < 0.

A frequency vector `w` can be optionally provided.

If `hz=true`, the plot x-axis will be displayed in Hertz, the input frequency vector is still treated as rad/s.

`kwargs` is sent as argument to Plots.plot.

See [`passivity_index`](@ref) for additional details.
See also [`ispassive`](@ref), [`passivity_index`](@ref).
"""
passivityplot
@recipe function passivityplot(p::Passivityplot; hz=false)
    systems, w = ControlSystemsBase._processfreqplot(Val{:sigma}(), p.args...)
    ws = (hz ? 1/(2π) : 1) .* w
    ny, nu = size(systems[1])
    nw = length(w)
    xguide --> (hz ? "Frequency [Hz]" : "Frequency [rad/s]")
    yguide --> "Passivity index"
    for (si, s) in enumerate(systems)
        s.ny == s.nu || throw(ArgumentError("passivity_index only defined for square systems"))
        Is = ss(I(s.ny), s.timeevol)
        G = (Is-s)feedback(Is, s)
        sv = sigma(G, w)[1] |> permutedims
        @series begin
            xscale --> :log10
            yscale --> :log10
            label --> "System $si"
            to1series(ws, sv)
        end
    end
    @series begin
        seriestype := :hline
        primary --> false
        linestyle --> :dash
        linecolor --> :black
        [1]
    end
end

"This is a helper function to make multiple series into one series separated by `Inf`. This makes plotting vastly more efficient."
function to1series(x::AbstractVector, y)
    r,c = size(y)
    y2 = vec([y; fill(Inf, 1, c)])
    x2 = repeat([x; Inf], c)
    x2,y2
end
