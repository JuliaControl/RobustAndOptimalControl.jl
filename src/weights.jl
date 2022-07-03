"""
    makeweight(low, f_mid, high)
    makeweight(low, (f_mid, gain_mid), high)

Create a weighting function that goes from gain `low` at zero frequency, through gain `gain_mid` to gain `high` at ∞

# Arguments:
- `low`: A number specifying the DC gain 
- `mid`: A number specifying the frequency at which the gain is 1, or a tuple `(freq, gain)`. If `gain_mid` is not specified, the geometric mean of `high` and `low` is used.
- `high`: A number specifying the gain at ∞

```@example
using ControlSystems, Plots
W = makeweight(10, (5,2), 1/10)
bodeplot(W)
hline!([10, 2, 1/10], l=(:black, :dash), primary=false)
vline!([5], l=(:black, :dash), primary=false)
```
"""
function makeweight(low, mid::Number, high)
    makeweight(low, (mid, high < 1 ? √(high*low) : 1), high)
end


function makeweight(low, mid, high)
    freq, mag = mid
    p = freq * √(((high/mag)^2-1)/(1-(low/mag)^2))
    W = tf([high, low*p], [1, p])
    W = ss(W)
end



"""
    neglected_delay(Lmax)

Return a multiplicative weight to represent the uncertainty coming from neglecting the dynamics `exp(-s*L)`
where `L ≤ Lmax`.
"Multivariable Feedback Control: Analysis and Design" Ch 7.4.5

See also [`gain_and_delay_uncertainty`](@ref) and [`neglected_lag`](@ref).

# Example:
```julia
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0) # Plant
W0 = neglected_delay(0.005) |> ss # Weight
W = I(2) + W0*I(2) * uss([δc(), δc()]) # Create a diagonal real uncertainty weighted in frequency by W0
Ps = P*W # Uncertain plant
Psamples = rand(Ps, 500) # Sample the uncertain plant for plotting
w = exp10.(LinRange(-1, 3, 300)) # Frequency vector
bodeplot(Psamples, w)
```
"""
function neglected_delay(Lmax)
    gain_and_delay_uncertainty(1, 1, Lmax)
end

"""
    gain_and_delay_uncertainty(kmin, kmax, Lmax)

Return a multiplicative weight to represent the uncertainty coming from neglecting the dynamics `k*exp(-s*L)`
where `k ∈ [kmin, kmax]` and `L ≤ Lmax`.
This weight is slightly optimistic, an expression for a more exact weight appears in eq (7.27), "Multivariable Feedback Control: Analysis and Design"

See also [`neglected_lag`](@ref) and [`neglected_delay`](@ref).

# Example:
```julia
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0) # Plant
W0 = gain_and_delay_uncertainty(0.5, 2, 0.005) |> ss # Weight
W = I(2) + W0*I(2) * uss([δc(), δc()]) # Create a diagonal real uncertainty weighted in frequency by W0
Ps = P*W # Uncertain plant
Psamples = rand(Ps, 500) # Sample the uncertain plant for plotting
w = exp10.(LinRange(-1, 3, 300)) # Frequency vector
bodeplot(Psamples, w)
```
"""
function gain_and_delay_uncertainty(kmin, kmax, Lmax)
    kb = (kmin+kmax)/2
    rk = (kmax-kmin)/2/kb
    tf([(1+rk/2)*Lmax, rk], [Lmax/2, 1])
end

"""
    neglected_lag(τmax)

Return a multiplicative weight to represent the uncertainty coming from neglecting the dynamics `1/(s*τ + 1)`
where `τ ≤ τmax`.
"Multivariable Feedback Control: Analysis and Design" Ch 7.4.5

See also [`gain_and_delay_uncertainty`](@ref) and [`neglected_delay`](@ref).

# Example:
```julia
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0) # Plant
W0 = neglected_lag(0.05) |> ss # Weight
W = I(2) + W0*I(2) * uss([δc(), δc()]) # Create a diagonal real uncertainty weighted in frequency by W0
Ps = P*W # Uncertain plant
Psamples = rand(Ps, 100) # Sample the uncertain plant for plotting
w = exp10.(LinRange(-1, 3, 300)) # Frequency vector
sigmaplot(Psamples, w)
```
"""
function neglected_lag(τmax)
    tf([τmax, 0], [τmax, 1])
end



"""
    centers, radii = fit_complex_perturbations(P, w; relative=true, nominal=:mean)

For each frequency in `w`, fit a circle in the complex plane that contains all models in the model set `P`, represented as an `LTISystem` with `Particles` coefficients. Note, the resulting radii can be quite unstable if the number of particles is small, in particular if the particles are normally distributed instead of uniformly.

If `realtive = true`, circles encompassing `|(P - Pn)/Pn|` will be returned (multiplicative/relative uncertainty).
If `realtive = false`, circles encompassing `|P - Pn|` will be returned (additive uncertainty).

If `nominal = :mean`, the mean of `P` will be used as nominal model. If `nominal = :first`, the first particle will be used. See `MonteCarloMeasurements.with_nominal` to set the nominal value in the first particle. 

See also [`nyquistcircles`](@ref) to plot circles (only if relative=false).
"""
function fit_complex_perturbations(P, w; relative=true, nominal=:mean)
    ControlSystems.issiso(P) || throw(ArgumentError("This function only works for SISO systems."))
    r = freqresp(P, w)
    centers = Vector{Complex{eltype(w)}}(undef, length(w))
    radii = Vector{eltype(w)}(undef, length(w))
    for i in axes(r,1)
        ri = r[i,1,1] # NOTE: assumes SISO
        c = nominal === :mean ? pmean(ri) : MonteCarloMeasurements.nominal(ri)
        rad = abs(relative ? pmaximum((ri - c)/c) : pmaximum(ri - c))
        centers[i] = c
        radii[i] = rad
    end
    centers, radii
end


"""
    nyquistcircles(w, centers, radii)

Plot the nyquist curve with circles. It only makes sense to call this function if the circles represent additive uncertainty, i.e., were calculated with `relative=false`.

See also [`fit_complex_perturbations`](@ref)
"""
nyquistcircles
@userplot nyquistcircles
@recipe function nyquistcircles(h::nyquistcircles)
    w, centers, radii = h.args
    t = LinRange(0, 2pi, 100)
    re, im = cos.(t), sin.(t)
    @series begin
        lalbe --> "Center"
        centers
    end
    for i in eachindex(radii)
        @series begin
            alpha --> 0.3
            color --> 1
            fillcolor --> 1
            label --> ""
            re .* radii[i] .+ centers[i].re, im .* radii[i] .+ centers[i].im
        end
    end
end

