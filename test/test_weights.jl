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
# `low` and `high` straddle 1, so the default gain_mid is 1 (preserved back-compat).
@test abs(evalfr(w, 1im)[]) ≈ 1 atol=1e-6

# Regression: previously, when both `low > 1` and `high > 1`, the default gain_mid was
# hard-coded to 1 which violates the formula's `low < mag < high` precondition and
# produced complex/NaN poles. Now falls back to √(low*high), which lies between them.
w_both_gt_1 = makeweight(2, 1, 3)
@test all(isreal, denvec(tf(w_both_gt_1))[1])
@test all(isreal, numvec(tf(w_both_gt_1))[1])
@test dcgain(w_both_gt_1)[] ≈ 2
@test evalfr(w_both_gt_1, 10000im)[] ≈ 3 atol=1e-3
@test abs(evalfr(w_both_gt_1, 1im)[]) ≈ √(2*3) atol=1e-6

# Symmetric same-side case: both below 1. Default gain_mid = √(low*high).
w_both_lt_1 = makeweight(0.5, 1, 0.1)
@test all(isreal, denvec(tf(w_both_lt_1))[1])
@test dcgain(w_both_lt_1)[] ≈ 0.5
@test evalfr(w_both_lt_1, 10000im)[] ≈ 0.1 atol=1e-3
@test abs(evalfr(w_both_lt_1, 1im)[]) ≈ √(0.5*0.1) atol=1e-6


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