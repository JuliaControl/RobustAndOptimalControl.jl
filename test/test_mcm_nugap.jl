using ControlSystemsBase
using RobustAndOptimalControl
using MonteCarloMeasurements

ω = with_nominal(0.9 .. 1.1, 1)
ζ = with_nominal(0.5 ± 0.01, 0.5)
G = tf([ω^2], [1, 2*ζ*ω, ω^2]) |> ss

@test nominal(G) == ss(tf([1], [1, 2*0.5*1, 1]))

gap = nugap(G)
@test gap isa typeof(ω)

# plot(gap)

# unsafe_comparisons(true)
# w = exp10.(LinRange(-2, 2, 200))
# bodeplot(G)

Gr = nu_reduction(G, 0.05)
@test 10 < nparticles(Gr.A) < nparticles(G.A)

Gr = nu_reduction_recursive(G, 0.05; verbose=true)
@test 1 < nparticles(Gr.A) < nparticles(G.A)

if isinteractive()
    unsafe_comparisons(true)
    bodeplot(G); bodeplot!(Gr, c=2)
    # bodeplot(Gr, c=2)
end


