using RobustAndOptimalControl, ControlSystems


G = ssrand(1,1,3, proper=true)
Gd = add_low_frequency_disturbance(G)
@test Gd.nx == 4
@test rank(obsv(Gd)) == 4 
@test rank(ctrb(Gd)) == 3
@test any(isapprox(0, atol=eps()), pole(Gd))

Gd = add_low_frequency_disturbance(G, measurement=true)
Gd.C[end] == 1
@test Gd.nx == 4
@test rank(obsv(Gd)) == 4
@test rank(ctrb(Gd)) == 3
@test any(isapprox(0, atol=eps()), pole(Gd))

G = ssrand(1,1,3, proper=true)
Gd = add_resonant_disturbance(G, 1, 0, 3)
@test Gd.nx == 5
@test rank(obsv(Gd)) == 5
@test rank(ctrb(Gd)) == 3
@test any(isapprox(1, atol=eps()), imag.(pole(Gd)))

Gd = add_resonant_disturbance(G, 1, 0, 1, measurement=true)
@test Gd.nx == 5
@test rank(obsv(Gd)) == 5
@test rank(ctrb(Gd)) == 3
@test any(isapprox(1, atol=eps()), imag.(pole(Gd)))



##

G = ssrand(1,1,3, proper=true, Ts=1)
GD = ssrand(1,1,3, proper=false, Ts=1)


## Diff
Gd = add_output_differentiator(G)
Gd2 = [tf(1,1); tf([1, -1], [1], 1)]*tf(G)
@test tf(Gd) ≈ Gd2
# @test hinfnorm(Gd-Gd2)[1] < 1e-10 hinfnorm not robust

## Int
Gd = add_output_integrator(G)
Gd2 = [tf(1,1); tf(1, [1, -1], 1)]*G
@test Gd ≈ Gd2
@test sminreal(Gd[1,1]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss

Gd = add_input_integrator(G)
@test sminreal(Gd[1,1]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss
@test tf(sminreal(Gd[2,1])) == tf(1, [1,-1], 1)

G = ssrand(2,1,3, proper=true, Ts=1)
Gd = add_input_integrator(G)
@test sminreal(Gd[1:2,1]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss
@test tf(sminreal(Gd[3,1])) == tf(1, [1,-1], 1)

G = ssrand(1,2,3, proper=true, Ts=1)
Gd = add_input_integrator(G)
@test sminreal(Gd[1,1:2]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss
@test tf(sminreal(Gd[2,1])) == tf(1, [1,-1], 1)
@test tf(sminreal(Gd[2,2])) == tf(0, 1)

Gd = add_input_integrator(G, 2) # other input
@test sminreal(Gd[1,1:2]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss
@test tf(sminreal(Gd[2,2])) == tf(1, [1,-1], 1)
@test tf(sminreal(Gd[2,1])) == tf(0, 1)

# Both
Gd = add_output_differentiator(G)
Gd = add_output_integrator(Gd, 1)
Gd2 = [tf(1,1); tf([1, -1], [1], 1); tf(1, [1, -1], 1)]*tf(G)
w = exp10.(LinRange(-2, 2, 100))

# These are harder to compare, the freqresp method seems most reliable
f1 = freqresp(Gd, w)
f2 = freqresp(Gd2, w)
@test f1 ≈ f2
# bodeplot([Gd, Gd2])