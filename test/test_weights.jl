using RobustAndOptimalControl
using ControlSystemsBase
using MonteCarloMeasurements

w = neglected_delay(1)
@test dcgain(w)[] < 1e-6
@test evalfr(w, 10000im)[] ≈ 2 atol=1e-3

w = gain_and_delay_uncertainty(1, 2, 1)


w = makeweight(0.1, 1, 2)
@test dcgain(w)[] ≈ 0.1
@test evalfr(w, 10000im)[] ≈ 2 atol=1e-3


w = neglected_lag(1)
@test dcgain(w)[] < 1e-6
@test evalfr(w, 10000im)[] ≈ 1 atol=1e-3


P = tf(1 ± 0.1, [1, 1 ± 0.2, 1])
w = 2π .* exp10.(LinRange(-1, 1, 200))
centers, radii = fit_complex_perturbations(P, w; relative=true, nominal=:mean)

@test 19 <= argmax(radii) <= 22

@test maximum(radii) > 1.4 
@test minimum(radii) < 0.5


centers, radii = fit_complex_perturbations(P, w; relative=true, nominal=:center)

@test maximum(radii) > 0.6
@test minimum(radii) < 0.5

##
w = 2π .* exp10.(LinRange(-2, 1, 200))
nyquistplot(P, w, ylims=(-2,2), xlims=(-2,2))

centers, radii = fit_complex_perturbations(P, w; relative=false, nominal=:mean)
nyquistcircles!(w, centers, radii)

centers, radii = fit_complex_perturbations(P, w; relative=false, nominal=:center)
nyquistcircles!(w, centers, radii, linecolor=:red)