# Uncertainty modeling
We provide two general means of modeling uncertainty, the traditional $M\Delta$ framework [^Skogestad][^Doyle91], and using parametric uncertainty. Support for parametric uncertainty is almost universal in Julia, not only in ControlSystems.jl, by means of computing with uncertain number types. In this tutorial, we will use a Monte-Carlo approach using uncertain number types from [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl).

Both the $M\Delta$ framework and parametric-uncertainty approaches are illustrated below.


```@contents
Pages = ["uncertainty.md"]
Depth = 3
```

## Uncertainty API
- [`δc`](@ref) Creates an uncertain complex parameter.
- [`δr`](@ref) Creates an uncertain real parameter.
- [`δss`](@ref) (Experimental) Creates an uncertain statespace model.
- [`neglected_delay`](@ref) Create a multiplicative weight that represents uncertainty from an unmodeled delay.
- [`neglected_lag`](@ref) Create a multiplicative weight that represents uncertainty from an unmodeled lag (pole).
- [`gain_and_delay_uncertainty`](@ref) Create a multiplicative weight that represents uncertainty from uncertain gains and delay.
- [`makeweight`](@ref) Create a custom weighting function.
- [`fit_complex_perturbations`](@ref)
- See [MonteCarloMeasurements.jl](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/) to create uncertain parameters that are represented by samples.

See example [`uncertain.jl`](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/uncertain.jl).


## Parametric uncertainty using [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl)
The most straightforward way to model uncertainty is to use uncertain parameters, using tools such as [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl) (strict, worst case guarantees) or [MonteCarloMeasurements](https://github.com/baggepinnen/MonteCarloMeasurements.jl) (less strict worst-case analysis or probabilistic).
In the following, we show an example with MIMO systems with both parametric uncertainty and diagonal, complex uncertainty, adapted from 8.11.3 in Skogestad, "Multivariable Feedback Control: Analysis and Design". This example is also available as a julia script in [`uncertain.jl`](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/uncertain.jl).

We will create uncertain parameters using the [`δr`](@ref) constructor from this package. One may alternatively create uncertain parameters directly using any of the constructors from MonteCarloMeasurements.jl. Most functions from ControlSystemsBase.jl should work with systems containing parameters from MonteCarloMeasurements.jl.

### Basic example
This example shows how to use MonteCarloMeasurements directly to build uncertain systems.
```@example BASIC_MCM
using ControlSystemsBase, MonteCarloMeasurements, Plots
ω = 1 ± 0.1 # Create an uncertain Gaussian parameter
```

```@example BASIC_MCM
ζ = 0.3..0.4 # Create an uncertain uniform parameter
```

```@example BASIC_MCM
G = tf(ω^2, [1, 2ζ*ω, ω^2]) # systems accept uncertain parameters
```

```@example BASIC_MCM
w = exp10.(-2:0.02:2)
bodeplot(G, w)
```

```@example BASIC_MCM
plot(step(G, 0:0.1:20))
```
### Example: Spinning satellite
This example makes use of real-valued uncertain parameters created using [`δr`](@ref), it comes from section 3.7.1 of Skogestad's book.
```@example satellite
using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements, Plots, LinearAlgebra
default(size=(640,480))
unsafe_comparisons(true)

a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)
K = ss(1.0I(2))

w = 2π .* exp10.(LinRange(-2, 2, 500))
S, PS, CS, T = gangoffour(P, K)
sigmaplot(S, w, lab="S")
sigmaplot!(T, w, c=2, lab="T", ylims=(0.01, 45))
```

Both sensitivity functions are very large, expect a non-robust system!

Next, we add parametric uncertainty
```@example satellite
a = 10*(1 + 0.1δr(100)) # Create an uncertain parameter with nominal value 10 and 10% uncertainty, represented by 100 samples
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)

Sp, PSp, CSp, Tp = gangoffour(P, K)
sigmaplot(Sp, w, lab="S")
sigmaplot!(Tp, w, c=2, lab="T", ylims=(0.01, 100))
```

Not only are sensitivity functions large, they vary a lot under the considered uncertainty. We can also plot a step response of one of the sensitivity functions to check how the system behaves
```@example satellite
plot(step(c2d(Tp, 0.01), 10))
```
This kind of plot is quite useful, it immediately tells you that this transfer function appears stable, and that there is uncertainty in the static gain etc.

Next, we add complex diagonal multiplicative input uncertainty. With input uncertainty of magnitude
$ϵ < \dfrac{1}{σ̄(T)}$ we are guaranteed robust stability (even for “full-block complex perturbations")

```@example satellite
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)

W0 = makeweight(0.2, (20,1), 2)
W = I(2) + W0 .* diagm([δc(100), δc(100)]) # Create a diagonal complex uncertainty weighted in frequency by W0, use 100 samples
Ps = P*W
Ss, PSs, CSs, Ts = gangoffour(Ps, K)
sigmaplot(Ss, w, lab="S")
sigmaplot!(Ts, w, c=2, lab="T", ylims=(0.01, 100))
```

Under this uncertainty, the sensitivity could potentially be sky high., note how some of the 100 realizations peak much higher than the others. This is an indication that the system might be unstable.

With complex entries in the system model, we can't really plot the step response, but we can plot, e.g., the absolute value
```@example satellite
res = step(c2d(Ts, 0.01), 10)
plot(res.t, [abs.(res.y)[1,:,1] abs.(res.y)[2,:,2]]) # plot only the diagonal response
```
Looks unstable to me. The analysis using $M\Delta$ methodology below will also reach this conclusion.


### Example: Distillation Process
This example comes from section 3.7.2 of Skogestad's book. In this example, we'll explore also complex uncertainties, created using [`δc`](@ref).
```@example distill
using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements, Plots, LinearAlgebra
default(size=(640,480))
unsafe_comparisons(true)

M = [87.8 -86.4; 108.2 -109.6]
G = Ref(ss(tf(1, [75, 1]))) .* M
RGA = relative_gain_array(G, 0)
sum(abs, RGA) # A good estimate of the true condition number, which is 141.7
```
large elements in the RGA indicate a process that is difficult to control

We consider the following inverse-based controller, which may also be looked upon as a steady-state decoupler with a PI controller
```@example distill
k1 = 0.7
Kinv = Ref(ss(tf(k1*[75, 1], [1, 0]))) .* inv(M) 

# reference filter
F = tf(1, [5, 1]) .* I(2)

w = 2π .* exp10.(LinRange(-2, 2, 500))
sigmaplot(input_sensitivity(G, Kinv), w)
sigmaplot!(output_sensitivity(G, Kinv), w, c=2)
```

Sensitivity looks nice, how about step response
```@example distill
plot(step(feedback(G*Kinv)*F, 20))
```

Looks excellent..

We consider again the input gain uncertainty as in the previous example, and we manually select the perturbations to be $ϵ_1 = 0.2$ and $ϵ_2 = 0.2$. We then have
```@example distill
G′ = G * diagm([1 + 0.2, 1 - 0.2])
plot!(step(feedback(G′*Kinv)*F, 20), l=:dash)
```

Looks very poor! The system was not robust to simultaneous input uncertainty!

We can also do this with a real, diagonal input uncertainty that grows with frequency
```@example distill
W0 = makeweight(0.2, 1, 2.0) # uncertainty goes from 20% at low frequencies to 200% at high frequencies
W = I(2) + W0 .* diagm([δr(100), δr(100)])
Gs = G*W

plot(step(feedback(G*Kinv)*F, 20))
plot!(step(feedback(G′*Kinv)*F, 20), l=:dash)
res = step(c2d(feedback(Gs*Kinv)*F, 0.01), 20)
mcplot!(res.t, abs.(res.y[:, :, 1]'), alpha=0.3)
mcplot!(res.t, abs.(res.y[:, :, 2]'), alpha=0.3)
```

The system is very sensitive to real input uncertainty!

With a complex, diagonal uncertainty, modeling both gain and phase variations, it looks slightly worse, but not much worse than with real uncertainty.
```@example distill
W = I(2) + W0 .* diagm([δc(100), δc(100)]) # note δc instead of δr above
Gs = G*W
res = step(c2d(feedback(Gs*Kinv)*F, 0.01), 20)
mcplot!(res.t, abs.(res.y[:, :, 1]'), alpha=0.3)
mcplot!(res.t, abs.(res.y[:, :, 2]'), alpha=0.3)
```

How about the sensitivity functions?
```@example distill
Si = input_sensitivity(Gs, Kinv)
sigmaplot(Si, w, c=1, lab="Si")
So = output_sensitivity(Gs, Kinv)
sigmaplot!(So, w, c=2, lab="So")
```

The sensitivity at the plant output is enormous. A low sensitivity with the nominal system does not guarantee robustness!

## Model-order reduction for uncertain models
The ``\nu``-gap metric is a measure of distance between models when they are used in a feedback loop. This metric has the nice property that a controller designed for a process ``P`` that achieves a normalized coprime factor margin ([`ncfmargin`](@ref)) of ``m``, will stabilize all models that are within a ``\nu``-gap distance of ``m`` from ``P``. This can be used to reduce the number of uncertain realizations for a model represented with `Particles` like above. Say that we have a plant model ``P``
```@example MCM_NUGAP
using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements
ω = with_nominal(0.9 .. 1.1, 1)
ζ = with_nominal(0.5 ± 0.01, 0.5)
P = tf([ω^2], [1, 2*ζ*ω, ω^2]) |> ss
```
represented by 2000 samples. We can compute the ``\nu``-gap metric between each realization in ``P`` and the nominal value (encoded using `with_nominal` above):
```@example MCM_NUGAP
gap = nugap(P)
```
The worst-case gap is:
```@example MCM_NUGAP
pmaximum(gap) # p for "particle" maximum
```
That means that if we design a controller for the nominal ``P`` without any uncertainty, and make sure that it achieves an [`ncfmargin`](@ref) of at least this value, it will stabilize all realizations in ``P``.

We can also reduce the number of realizations in ``P`` by discarding those that are close in the ``\nu``-gap sense to the nominal value:
```@example MCM_NUGAP
Pr = nu_reduction(P, 0.1)
```
here, all realizations that were within a ``\nu``-gap distance of 0.1 from the nominal value were discarded. [`nu_reduction`](@ref) usually reduces the number of realizations substantially.

We can reduce the number of realizations even further using [`nu_reduction_recursive`](@ref):
```@example MCM_NUGAP
Gr = nu_reduction_recursive(P, 0.1)
```
The algorithm used in [`nu_reduction_recursive`](@ref) has a worst-case complexity of ``O(N^2)`` where ``N`` is the number of realizations (particles) in ``P``, but this complexity is only problematic for small gaps and large number of realizations, say, ``\nu < 0.02`` and ``N > 50``.

!!! warn "Stochastic interpretation"
     If `P` has a stochastic interpretation, i.e., the coefficients come from some distribution, this interpretation will be lost after reduction, mean values and standard deviations will not be preserved. The reduced system should instead be interpreted as preserving worst-case uncertainty.

The reduced model is useful for plotting and manual loop-shaping etc.

## Using the $M\Delta$ framework
The examples above never bothered with things like the "structured singular value", $\mu$ or linear-fractional transforms. We do, however, provide some elementary support for this modeling framework.



In robust control, we often find ourselves having to consider the feedback interconnections below.
```
        ┌─────────┐
  zΔ◄───┤         │◄────wΔ
        │         │
   z◄───┤    P    │◄────w
        │         │
   y◄───┤         │◄────u
        └─────────┘
```
```
        ┌─────────┐
  zΔ◄───┤         │◄────wΔ
        │         │
   z◄───┤    P    │◄────w
        │         │
   y┌───┤         │◄───┐u
    │   └─────────┘    │
    │      ┌───┐       │
    └─────►│ K ├───────┘
           └───┘
```
```
           ┌───┐
    ┌─────►│ Δ ├───────┐
    │      └───┘       │
    │   ┌─────────┐    │
  zΔ└───┤         │◄───┘wΔ
        │         │
   z◄───┤    P    │◄────w
        │         │
   y┌───┤         │◄───┐u
    │   └─────────┘    │
    │      ┌───┐       │
    └─────►│ K ├───────┘
           └───┘
```

The first block diagram denotes an open-loop system $P$ with an uncertainty mapping $w_\Delta = \Delta  z_\Delta$, a *performance mapping* from $w$ to $z$ and a input-output mapping between $u$ and $y$. Such a system $P$ can be partitioned as
```math
P = \begin{bmatrix}
P_{11} & P_{12} & P_{13}\\
P_{21} & P_{22} & P_{23}\\
P_{31} & P_{32} & P_{33}\\
\end{bmatrix}
```
where each $P(s)_{ij}$ is a transfer matrix. The type [`UncertainSS`](@ref) with constructor [`uss`](@ref) represents the block
```math
P = \begin{bmatrix}
P_{11} & P_{12}\\
P_{21} & P_{22}\\
\end{bmatrix}
```
while an [`ExtendedStateSpace`](@ref) object represents the block
```math
P = \begin{bmatrix}
P_{22} & P_{23}\\
P_{32} & P_{33}\\
\end{bmatrix}
```
there is thus no type that represents the full system $P$ above. However, we provide the function [`partition`](@ref) which allows you to convert from a regular statespace system to an extended statespace object, and it is thus possible to represent $P$ by placing the whole block 
```math
P = \begin{bmatrix}
P_{22} & P_{23}\\
P_{32} & P_{33}\\
\end{bmatrix}
```
into $P_{22}$ for the purposes of uncertainty analysis (use `ss` to convert it to a standard statespace object), and later use [`partition`](@ref) to recover the internal block structure. 

Given an [`UncertainSS`](@ref) $P$, we can close the loop around $\Delta$ by calling `lft(P, Δ, :u)`, and given an [`ExtendedStateSpace`](@ref), we can close the loop around `K` by calling `starprod(P, K)` or `lft(P, K)` (using positive feedback). This works even if `P` is a regular statespace object, in which case the convention is that the inputs and outputs are ordered as in the block diagrams above. The number of signals that will be connected by [`lft`](@ref) is determined by the input-output arity of $K$ and $\Delta$ respectively.

We have the following methods for `lft` (in addition to the standard ones in ControlSystemsBase.jl)
- `lft(G::UncertainSS, K::LTISystem)` forms the lower LFT closing the loop around $K$.
- `lft(G::UncertainSS, Δ::AbstractArray=G.Δ)` forms the upper LFT closing the loop around $\Delta$.
- `lft(G::ExtendedStateSpace, K)` forms the lower LFT closing the loop around $K$.

### Robust stability and performance
To check robust stability of the system in the last block diagram (with or without $z$ and $w$), we can use the functions [`structured_singular_value`](@ref), [`robstab`](@ref) and [`diskmargin`](@ref).

Currently, [`structured_singular_value`](@ref) is rather limited and supports diagonal complex blocks only. If $\Delta$ is a single full complex block, `opnorm(P.M) < 1` is the condition for stability.

Robust performance can be verified by introducing an additional fictitious "performance perturbation" $\Delta_p$ which is a full complex block, around which we close the loop from $z$ to $w$ and check the [`structured_singular_value`](@ref) with the augmented perturbation block
```math
\Delta_a = \begin{bmatrix}
\Delta & 0\\
0      & \Delta_p\\
\end{bmatrix}
```



### Examples
We repeat the first example here, but using $M\Delta$ formalism rather than direct Monte-Carlo modeling.

When we call [`δc`](@ref) without any arguments, we get a symbolic (or structured) representation of the uncertainty rather than the sampled representation we got from calling `δc(100)`.


```@example satellite
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)
W0 = makeweight(0.2, (1,1), 2) |> ss
W = I(2) + W0 .* uss([δc(), δc()]) # Create a diagonal complex uncertainty weighted in frequency by W0
Ps = P*W
```
`Ps` is now represented as a upper linear fractional transform (upper LFT).

We can draw samples from this uncertainty representation (sampling of $\Delta$ and closing the loop `starprod(Δ, Ps)`) like so
```@example satellite
Psamples = rand(Ps, 100)
sigmaplot(Psamples, w)
```

We can extract the nominal model using

```@example satellite
system_mapping(Ps)
```
And obtain $M$ and $\Delta$ when the loop is closed with $K$ has
```@example satellite
lft(Ps, K).M
```
```@example satellite
Ps.Δ # Ps.delta also works
```
We can evaluate the frequency response of $M$ and calculate the structured singular value $\mu$

```@example satellite
M = freqresp(lft(Ps, -K).M, w) # -K to get negative feedback
μ = structured_singular_value(M)
plot(w, μ, xscale=:log10)
```

$\mu$ is very high, whenever $\mu > 1$, the system is not stable with respect to the modeled uncertainty.
The tolerated uncertainty is only about $\dfrac{1}{||\mu||_\infty}$
```@example satellite
1/norm(μ, Inf)
```
of the modeled uncertainty. Another way of calculating this value is
```@example satellite
robstab(lft(Ps, -K))
```




### Internals of the $M\Delta$ framework
TODO


[^Skogestad]: Skogestad, "Multivariable Feedback Control: Analysis and Design"

[^Doyle91]: Doyle, Packard, Zhou, "Review of LFTs, LMIs and μ". [`https://www.researchgate.net/publication/257200344_Review_of_LFTs_LMIs_and_mu`](https://www.researchgate.net/publication/257200344_Review_of_LFTs_LMIs_and_mu)


## Uncertain time delays

Modeling uncertain time delays can be done in several ways, one approach is to make use of a multiplicative uncertainty weight created using [`neglected_delay`](@ref) multiplied by an uncertain element created using [`δc`](@ref), example:
```@example uncertain_delay
using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements, Plots, LinearAlgebra
a  = 10
P  = ss([0 a; -a 0], I(2), [1 a; -a 1], 0) # Plant
W0 = neglected_delay(0.005) |> ss # Weight
W  = I(2) + W0 .* uss([δc(), δc()]) # Create a diagonal real uncertainty weighted in frequency by W0
Ps = P*W # Uncertain plant
Psamples = rand(Ps, 500) # Sample the uncertain plant for plotting
w = exp10.(LinRange(-1, 3, 300)) # Frequency vector
bodeplot(Psamples, w, legend=false)
```
Note how this approximation approach imparts some uncertainty also in the gain.

More details on this approach can be found in Skogestad sec. 7.4.

The other alternative is to use use sampled uncertain delays. The next example shows how we can create a system with an uncertain delay, where we know that the delay is an integer number of milliseconds between 1ms and 4ms.
```@example uncertain_delay
using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements, Plots, LinearAlgebra
unsafe_comparisons(true)
L = Particles(collect((1:4) ./ 1000)) # Uncertain time delay, an integer number of milliseconds between 1ms and 4ms
P = delay(L)*tf(1, [0.01, 1])
C = pid(2, 1)
w = exp10.(-1:0.01:4)
plot(
     bodeplot(P, exp10.(-1:0.001:3), legend=false),
     # plot(step(feedback(P, C), 0:0.0001:0.05), lab="L = " .* string.(P.Tau[].particles'), title="Disturbance response"), # This simulation requires using ControlSystems
     nyquistplot(P*C, w[1:10:end], points=true, xlims=(-3.5, 2.5), ylims=(-5, 1.5), Ms_circles=[1.5, 2], alpha=1) # Note, the nyquistplot with uncertain coefficients requires manual selection of plot limits
)
```
Notice how the gain is completely certain, while the phase starts becoming very uncertain for high frequencies.