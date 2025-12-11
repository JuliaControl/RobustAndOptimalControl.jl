using RobustAndOptimalControl, ControlSystemsBase


G = ssrand(1,1,3, proper=true)
Gd = add_low_frequency_disturbance(G)
@test Gd.nx == 4
@test rank(obsv(Gd)) == 4 
@test rank(ctrb(Gd)) == 3
@test any(isapprox(0, atol=eps()), poles(Gd))

Gd = add_low_frequency_disturbance(G, measurement=true)
Gd.C[end] == 1
@test Gd.nx == 4
@test rank(obsv(Gd)) == 4
@test rank(ctrb(Gd)) == 3
@test any(isapprox(0, atol=eps()), poles(Gd))

G = ssrand(4,2,3, proper=true)
Gd = add_low_frequency_disturbance(G)
@test Gd.nx == G.nx + G.nu
@test rank(obsv(Gd)) == Gd.nx
@test rank(ctrb(Gd)) == G.nx
@test any(isapprox(0, atol=eps()), poles(Gd))
@test Gd.A[end-1:end, end-1:end] == 0I

G = ssrand(2,4,3, proper=true)
Gd = add_low_frequency_disturbance(G, measurement=true)
@test Gd.nx == G.nx + G.ny
@test rank(obsv(Gd)) == G.nx + G.ny
@test rank(ctrb(Gd)) == G.nx
@test any(isapprox(0, atol=eps()), poles(Gd))
@test Gd.C[:, end-1:end] == I
@test Gd.A[end-1:end, end-1:end] == 0I

# Discrete time
G = ssrand(4,2,3, proper=true, Ts=0.1)
Gd = add_low_frequency_disturbance(G)
@test Gd.nx == G.nx + G.nu
@test rank(obsv(Gd)) == Gd.nx
@test rank(ctrb(Gd)) == G.nx
@test any(isapprox(1, atol=eps()), poles(Gd))
@test Gd.A[end-1:end, end-1:end] == I

G = ssrand(2,4,3, proper=true, Ts=0.1)
Gd = add_low_frequency_disturbance(G, measurement=true)
@test Gd.nx == G.nx + G.ny
@test rank(obsv(Gd)) == G.nx + G.ny
@test rank(ctrb(Gd)) == G.nx
@test any(isapprox(1, atol=eps()), poles(Gd))
@test Gd.C[:, end-1:end] == I
@test Gd.A[end-1:end, end-1:end] == I


G = ssrand(1,1,3, proper=true)
Gd = add_resonant_disturbance(G, 1, 0, 3)
@test Gd.nx == 5
@test rank(obsv(Gd)) == 5
@test rank(ctrb(Gd)) == 3
@test any(isapprox(1, atol=eps()), imag.(poles(Gd)))

Gd = add_resonant_disturbance(G, 1, 0, 1, measurement=true)
@test Gd.nx == 5
@test rank(obsv(Gd)) == 5
@test rank(ctrb(Gd)) == 3
@test any(isapprox(1, atol=eps()), imag.(poles(Gd)))

# Discrete time and input matrix
G = c2d(ss(tf(1.0, [1, 1])), 0.1)
Gd = add_resonant_disturbance(G, 1, 0, [1.0])
@test sminreal(Gd) == G
@test Gd.nx == 3
allapproxin(a, b) = all(any(a .≈ b', dims=2))
@test allapproxin(poles(Gd), [eigvals(exp([0 -1; 1 0]*0.1)); exp(-1*0.1)])

# Discrete time and input index
G = c2d(ss(tf(1.0, [1, 1])), 0.1)
Gd = add_resonant_disturbance(G, 1, 0, 1)
@test sminreal(Gd) == G
@test Gd.nx == 3
allapproxin(a, b) = all(any(a .≈ b', dims=2))
@test allapproxin(poles(Gd), [eigvals(exp([0 -1; 1 0]*0.1)); exp(-1*0.1)])

# Discrete time, input matrix, and measurement disturbance
G = c2d(ss(tf(1.0, [1, 1])), 0.1)
Gd = add_resonant_disturbance(G, 1, 0, [1.0 0.0], measurement=true)
@test sminreal(Gd) == G
@test Gd.nx == 3
allapproxin(a, b) = all(any(a .≈ b', dims=2))
@test allapproxin(poles(Gd), [eigvals(exp([0 -1; 1 0]*0.1)); exp(-1*0.1)])
@test size(Gd.C, 2) == 3  # C matrix extended with disturbance states


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
@test tf(Gd) ≈ tf(Gd2)
@test sminreal(Gd[1,1]) == G # Exact equivalence should hold here
@test Gd.nx == 4 # To guard agains changes in realization of tf as ss


Gc = ssrand(1,1,3, proper=true)
Gdc = add_output_integrator(Gc)
Gd2c = [tf(1); tf(1, [1, 0])]*Gc
@test tf(Gdc) ≈ tf(Gd2c)
@test sminreal(Gdc[1,1]) == Gc # Exact equivalence should hold here
@test Gdc.nx == 4 # To guard agains changes in realization of tf as ss

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


# Input diff
nx,nu,ny = G.nx, G.nu, G.ny
Gd = add_input_differentiator(G, 1:G.nu)
@test Gd.A[nx+1:end, nx+1:end] == 0I
@test Gd.B[nx+1:end, :] == I
@test Gd.C[ny+1:end, nx+1:end] == -I
@test Gd.D[ny+1:end, :] == I