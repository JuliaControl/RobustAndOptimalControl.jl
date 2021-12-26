using RobustAndOptimalControl, ControlSystems, MonteCarloMeasurements, Plots
unsafe_comparisons(true)

using RobustAndOptimalControl: δr, δc

# example in section 3.7.1, spinning satellite
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)
K = ss(I(2))

w = 2π .* exp10.(LinRange(-2, 2, 500))
S, PS, CS, T = RobustAndOptimalControl.gangoffour2(P, K)
sigmaplot(S, w, lab="S")
sigmaplot!(T, w, c=2, lab="T", ylims=(0.01, 45))
# Both sensitivity functions are very large, expect a non-robust system!

## Parametric uncertainty
a = 10*(1 + δr(100))
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)

Sp, PSp, CSp, Tp = RobustAndOptimalControl.gangoffour2(P, K)
sigmaplot(Sp, w, lab="S")
sigmaplot!(Tp, w, c=2, lab="T", ylims=(0.01, 45))
# Not only are sensitivity functions large, they vary a lot under the considered uncertainty



## Add complex diagonal multiplicative input uncertainty
# With input uncertainty of magnitude ϵ < 1 / σ̄(T) we are guaranteed robust stability (even for “full-block complex perturbations")
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)

W0 = Weights.makeweight(0.2, 2)
W = I(2) + W0 * diagm([δc(100), δc(100)])
Ps = P*W
Ss, PSs, CSs, Ts = RobustAndOptimalControl.gangoffour2(Ps, K)
sigmaplot(Ss, w, lab="S")
sigmaplot!(Ts, w, c=2, lab="T", ylims=(0.01, 45))
# Under this uncertainty, the sensitivity could potentially be sky high., note how some of the 100 realizations peak much higher than the others



##

# 3.7.2 Motivating robustness example no. 2: Distillation Process
M = [87.8 -86.4; 108.2 -109.6]
G = ss(tf(1, [75, 1])) * M
RGA = relative_gain_array(G, 0)
sum(abs, RGA) # A good estimate of the true condition number, which is 141.7
# large elments in the RGA indicate a process that is difficult to control

# We consider the following inverse-based controller, which may also be looked
# upon as a steady-state decoupler with a PI controlle
k1 = 0.7
Kinv = ss(tf(k1*[75, 1], [1, 0])) * inv(M) 

# reference filter
F = tf(1, [5, 1])

sigmaplot(input_sensitivity(G, Kinv), w)
sigmaplot!(output_sensitivity(G, Kinv), w, c=2)
# Sensitivity looks nice, how about step response
plot(step(feedback(G*Kinv)*F, 20))
# Looks excellent..

# We consider again the input gain uncertainty as in the previous example,
# and we select ϵ1 = 0.2 and ϵ2 = 0.2. We then have
G′ = G * diagm([1 + 0.2, 1 - 0.2])
plot!(step(feedback(G′*Kinv)*F, 20), l=:dash)
# Looks very poor! The system was not robust to simultaneous input uncertainty!

# We can also do this with a complex, diagonal input uncertainty that grows with frequency
W0 = Weights.makeweight(0.2, 1, 2) # uncertainty goes from 20% at low frequencies to 200% at high frequencies
W = I(2) + W0 * diagm([δc(100), δc(100)])
Gs = G*W
res = step(c2d(feedback(Gs*Kinv)*F, 0.01), 20)
mcplot!(res.t, abs.(res.y[:, :, 1]'), alpha=0.3)
mcplot!(res.t, abs.(res.y[:, :, 2]'), alpha=0.3)
# The system is very sensitive to complex input uncertainty

# How about real input uncertainty?
W = I(2) + W0 * diagm([δr(100), δr(100)])
Gs = G*W

plot(step(feedback(G*Kinv)*F, 20))
plot!(step(feedback(G′*Kinv)*F, 20), l=:dash)
res = step(c2d(feedback(Gs*Kinv)*F, 0.01), 20)
mcplot!(res.t, abs.(res.y[:, :, 1]'), alpha=0.3)
mcplot!(res.t, abs.(res.y[:, :, 2]'), alpha=0.3)
# Also bad

Si = input_sensitivity(Gs, Kinv)
sigmaplot(Si, w, c=1, lab="Si")
So = output_sensitivity(Gs, Kinv)
sigmaplot!(So, w, c=2, lab="So")
# The sensitivity at the plant output is enormous. A low sensitivity with the nominal system does not guarantee robustness!