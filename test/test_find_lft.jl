using MonteCarloMeasurements
unsafe_comparisons()
u1 = 2*StaticParticles(32, Uniform(-1,1))
u2 = 5*StaticParticles(32, Uniform(-1,1))
u3 = 0.2*StaticParticles(32, Uniform(-1,1))
nx, nu, ny = 10, 5, 5
A = randn(nx,nx) .+ 0 .* StaticParticles.(32)
A[1,1] = 4 + u1
A[2, 3] = 9 + u2
B = randn(nx, nu) .+ 0 .* StaticParticles.(32)
B[4, 1] = 2 + u3

sys = ss(A, B, randn(ny,nx), 0)
l, res = find_lft(sys, [u1, u2, u3])

sysh = ss(l,sys)

@test mean.(sysh.A) ≈ mean.(sys.A) atol=1e-6
@test std.(sysh.A) ≈ std.(sys.A) atol=1e-6
@test mean.(sysh.B) ≈ mean.(sys.B) atol=1e-6
@test std.(sysh.B) ≈ std.(sys.B) atol=1e-6
@test mean.(sysh.C) ≈ mean.(sys.C) atol=1e-6
@test std.(sysh.C) ≈ std.(sys.C) atol=1e-6
@test mean.(sysh.D) ≈ mean.(sys.D) atol=1e-6
@test std.(sysh.D) ≈ std.(sys.D) atol=1e-6
