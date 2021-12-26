"""
    Diskmargin

The notation follows "An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet

# Fields:
`α`: The disk margin
`ω0`: The worst-case frequency
`f0`: The destabilizing perturbation `f0` is a complex number with simultaneous gain and phase variation. This critical perturbation causes an instability with closed-loop pole on the imaginary axis at the critical frequency ω0 
`δ0`: The uncertain element generating f0.
`γmin`: The lower real-axis intercept of the disk (classical lower gain margin).
`γmax`: The upper real-axis intercept of the disk (classical upper gain margin).
`ϕm`: is the classical phase margin.
`σ`: The skew parameter that was used to calculate the margin

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

function Base.show(io::IO, dm::Diskmargin)
    println(io, "Disk margin with:")
    println(io, "Margin: ", dm.α)
    println(io, "Frequency: ", dm.ω0)
    if dm.γmax < dm.γmin # In this case, we have an "inverted circle"
        println(io, "Gain margins: [$(dm.γmin), Inf]")
    else
        println(io, "Gain margins: [$(dm.γmin), $(dm.γmax)]")
    end
    println(io, "Phase margin: ", dm.ϕm)
    println(io, "Skew: ", dm.σ)
    println(io, "Worst-case perturbation: ", dm.f0)
end

"""
    Disk

Represents a perturbation disc in the complex plane. `Disk(0.5, 2)` represents all perturbations in the circle centered at 1.25 with radius 0.75, or in other words, a gain margin of 2 and a pahse margin of 36.9 degrees.

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
        ϕm = ϕm >= 1 ? Inf : rad2deg(acos(ϕm))
    end
    Disk(γmin, γmax, c, r, ϕm)
end

function Disk(; α, σ)
    γmin = (2 - α*(1-σ)) / (2 + α*(1+σ))
    γmax = (2 + α*(1-σ)) / (2 - α*(1+σ))
    Disk(γmin, γmax)
end

Disk(dm::Diskmargin) = Disk(dm.γmin, dm.γmax)
ControlSystems.nyquist(d::Disk) = Disk(-inv(d.γmin), -inv(d.γmax)) # translate the disk to a nyquist exclusion disk

"""
    diskmargin(L, σ = 0)
    diskmargin(L, σ::Real, ω)

Calculate the disk margin of LTI system `L`. `L` is supposed to be a loop-transfer function, i.e., it should be square. If `L = PC` the disk margin for output perturbations is computed, whereas if `L = CP`, input perturbations are considered. If the method `diskmargin(P, C, args...)` is used, both are computed.

The implementation and notation follows
"An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet
https://arxiv.org/abs/2003.04771

The margins are aviable as fields of the returned objects, see [`Diskmargin`](@ref).

# Arguments:
- `L`: A loop-transfer function.
- `σ`: If little is known about the distribution of gain variations then σ = 0
is a reasonable choice as it allows for a gain increase or decrease by the same relative amount.
The choice σ < 0 is justified if the gain can decrease by a larger factor than it can increase.
Similarly, the choice σ > 0 is justified when the gain can increase by a larger factor than it can
decrease.
If σ = −1 then the disk margin condition is αmax = inv(MT). This margin is related to the robust
stability condition for models with multiplicative uncertainty of the form P (1 + δ).
If σ = +1 then the disk margin condition is αmax = inv(MS)
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
plot(w, dms)
```
"""
function diskmargin(L::LTISystem, σ::Real=0; kwargs...)
    issiso(L) || return sim_diskmargin(L, σ)
    S̄ = 1/(1 + L) + (σ-1)/2
    n,ω0 = hinfnorm(S̄; kwargs...)
    diskmargin(L, σ, ω0)
end

diskmargin(L::LTISystem, σ::Real, ω::AbstractArray) = map(w->diskmargin(L, σ, w), ω)

function diskmargin(L::LTISystem, σ::Real, ω0::Real)
    issiso(L) || return sim_diskmargin(L, σ, [ω0])[]
    S̄ = 1/(1 + L) + (σ-1)/2
    freq = isdiscrete(L) ? cis(ω0*L.Ts) : complex(0, ω0)
    Sω = evalfr(S̄, freq)[]
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

function γϕcurve(α, σ; N = 200)
    θ = LinRange(0, π, N)
    f = @. (2 - α*cis(θ)*(1-σ)) / (2 + α*cis(θ)*(1+σ))
    @. (abs(f), abs(rad2deg(angle(f))))
end

@recipe function plot(d::Disk; nyquist = false)
    θ = LinRange(0, 2pi, 200)
    re, im = @. cos(θ), sin(θ)
    if nyquist
        d = ControlSystems.nyquist(d)
    end
    c,r = d.c, d.r
    @series begin 
        fill --> true
        fillalpha --> 0.5
        @. r*re+c, r*im
    end
end

@recipe function plot(dm::Diskmargin; nyquist=false)
    if nyquist
        @series begin
            label --> "σ = $(dm.σ)"
            ControlSystems.nyquist(Disk(dm))
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

@recipe function plot(dm::AbstractVector{<:Diskmargin})
    w = [dm.ω0 for dm in dm]
    layout --> (2, 1)
    link --> :x
    @series begin
        subplot --> 1
        title --> "Gain margin"
        label --> ["Lower" "Upper"]
        # xguide --> "Frequency"
        xscale --> :log10
        yscale --> :log10
        gma = getfield.(dm, :γmax)
        gmi = getfield.(dm, :γmin)
        # ylims --> (0, min(10, maximum(gma)))
        data = [gmi gma]
        data = max.(0, data)
        replace!(data, 0 => -Inf)
        w, data
    end
    @series begin
        subplot --> 2
        title --> "Phase margin"
        xguide --> "Frequency"
        xscale --> :log10
        label --> ""
        ϕm = getfield.(dm, :ϕm)
        w, ϕm
    end
end

@recipe function plot(dm::AbstractVector{<:AbstractVector{Diskmargin}})
    reduce(hcat, dm)
end


@recipe function plot(dm::AbstractMatrix{Diskmargin})
    w = getfield.(dm[:,1], :ω0)
    ny = size(dm, 2)
    length(w) == size(dm, 1) || throw(ArgumentError("Frequency vector and diskmargin vector must have the same lengths."))
    layout --> (2, 1)
    link --> :x
    neg2inf(x) = x <= 0 ? Inf : x
    for i = 1:ny
        is = ny == 1 ? "" : "$i"
        @series begin
            subplot --> 1
            title --> "Gain margin"
            # label --> permutedims([["Lower $i" for i = 1:ny]; ["Upper $i" for i = 1:ny]])
            label --> ["Upper"*is "Lower"*is]
            # xguide --> "Frequency"
            xscale --> :log10
            yscale --> :log10
            gma = getfield.(dm[:, i], :γmax) .|> neg2inf
            gmi = getfield.(dm[:, i], :γmin) .|> neg2inf
            replace!(gmi, 0 => -Inf)
            # ylims --> (0, min(10, maximum(gma)))
            w, [gmi gma]
        end
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