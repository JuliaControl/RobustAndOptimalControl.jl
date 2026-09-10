using RobustAndOptimalControl, ControlSystemsBase


G = ssrand(1,1,3, proper=true)
Gd = add_low_frequency_disturbance(G)
@test Gd.nx == 4
@test rank(obsv(Gd)) == 4 
@test rank(ctrb(Gd)) == 3
@test any(isapprox(0, atol=eps()), poles(Gd))

Gd = add_low_frequency_disturbance(G, measurement=true)
@test Gd.C[end] == 1
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

# Integer-Ai variant on a MIMO plant (nu > 1) — regression for sizing of Ad
G = ssrand(2, 3, 4, proper=true)
Gd = add_low_frequency_disturbance(G, 2)
@test Gd.nx == G.nx + 1
@test rank(obsv(Gd)) == Gd.nx
@test any(isapprox(0, atol=eps()), poles(Gd))
@test Gd.A[1:G.nx, end] == [0, 1, 0, 0]
@test Gd.A[end, end] == 0

# Same, discrete time
Gddisc = add_low_frequency_disturbance(c2d(G, 0.1), 2)
@test Gddisc.A[end, end] == 1

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

# Two-column Bd (both resonant states inject into the plant state)
G = c2d(ss(tf(1.0, [1, 1])), 0.1)
Gd = add_resonant_disturbance(G, 1, 0, [1.0 0.5])
@test Gd.nx == 3
@test allapproxin(poles(Gd), [eigvals(exp([0 -1; 1 0]*0.1)); exp(-1*0.1)])

# Input validation for Bd shape
@test_throws ArgumentError add_resonant_disturbance(G, 1, 0, zeros(G.nx, 3))
@test_throws ArgumentError add_resonant_disturbance(G, 1, 0, zeros(G.nx + 1, 1))
@test_throws ArgumentError add_resonant_disturbance(G, 1, 0, zeros(G.ny, 1), measurement=true)


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

# One integrator output and one integrator state per requested index, in the order given
Gm = ssrand(3,2,2, proper=true)
w = exp10.(LinRange(-2, 2, 100))
for inds in ([1], [2,3], [3,1], 1:3)
    Gi = add_output_integrator(Gm, inds)
    @test Gi.ny == Gm.ny + length(inds)
    @test Gi.nx == Gm.nx + length(inds)
    @test sminreal(Gi[1:Gm.ny, :]) == Gm # The original outputs are untouched
    @test Gi.A[Gm.nx+1:end, 1:Gm.nx] ≈ Gm.C[inds, :] # State k integrates output inds[k]
    # freqresp is the reliable comparison here, the tf of the augmented system carries an
    # uncancelled pole/zero pair at the origin for the non-integrated outputs
    @test freqresp(Gi[Gm.ny+1:end, :], w) ≈ freqresp(ss(tf(1, [1, 0])) .* Gm[inds, :], w)
end
# An integer index is equivalent to the length-one vector
@test add_output_integrator(Gm, 2) == add_output_integrator(Gm, [2])
@test_throws ArgumentError add_output_integrator(Gm, 4)
# `neg` negates the added outputs and leaves the integrator state dynamics untouched
Gi = add_output_integrator(Gm, [2,3])
Gin = add_output_integrator(Gm, [2,3]; neg=true)
@test Gin.A == Gi.A
@test Gin.B == Gi.B
@test Gin.C == [Gi.C[1:Gm.ny, :]; -Gi.C[Gm.ny+1:end, :]]

# The discrete integrator is a forward-Euler time integral, xᵢ⁺ = (1-ϵ)xᵢ + Ts*y
ϵ = 1e-3
Gmd = ssrand(2,2,2, proper=true, Ts=0.1)
Gid = add_output_integrator(Gmd, [2,1]; ϵ)
@test Gid.A[3:4, 1:2] ≈ Gmd.Ts * Gmd.C[[2,1], :]
@test Gid.A[3:4, 3:4] ≈ (1 - ϵ)*I(2)
wd = exp10.(LinRange(-2, 1, 100))
@test freqresp(Gid[3:4, :], wd) ≈ freqresp(ss(tf(Gmd.Ts, [1, -(1-ϵ)], Gmd.Ts)) .* Gmd[[2,1], :], wd)

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