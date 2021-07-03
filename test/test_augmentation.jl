using RobustAndOptimalControl, ControlSystems
G = ssrand(1,1,3, proper=true, Ts=1)
GD = ssrand(1,1,3, proper=false, Ts=1)


## Diff
Gd = add_differentiator(G)
Gd2 = [tf(1,1); tf([1, -1], [1], 1)]*tf(G)
tf(Gd) ≈ Gd2
@test hinfnorm(Gd-Gd2)[1] < 1e-10

## Int
Gd = add_integrator(G)
Gd2 = [tf(1,1); tf(1, [1, -1], 1)]*G
@test Gd ≈ Gd2
@test sminreal(Gd[1,1]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss

# Both
Gd = add_differentiator(G)
Gd = add_integrator(Gd, 1)
Gd2 = [tf(1,1); tf([1, -1], [1], 1); tf(1, [1, -1], 1)]*tf(G)
w = exp10.(LinRange(-2, 2, 100))

# These are harder to compare, the freqresp method seems most reliable
f1 = freqresp(Gd, w)
f2 = freqresp(Gd2, w)
@test f1 ≈ f2
# bodeplot([Gd, Gd2])