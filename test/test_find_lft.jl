using RobustAndOptimalControl, Random, Statistics, MonteCarloMeasurements, ControlSystems
Random.seed!(0)
unsafe_comparisons()
u1 = 2*δ()
u2 = 5*δ()
u3 = 0.2*δ()
nx, nu, ny = 10, 5, 5
A = randn(nx,nx) .+ 0 .* δ.()
A[1,1] = 4 + u1
A[2, 3] = 9 + u2
B = randn(nx, nu) .+ 0 .* δ.()
B[4, 1] = 2 + u3

sys = ss(A, B, randn(ny,nx), 0)
l, res = find_lft(sys, [u1, u2, u3])

sysh = ss(l,sys)

@test mean.(sysh.A) ≈ mean.(sys.A) atol=1e-5
@test std.(sysh.A) ≈ std.(sys.A) rtol=1e-5
@test mean.(sysh.B) ≈ mean.(sys.B) atol=1e-5
@test std.(sysh.B) ≈ std.(sys.B) rtol=1e-5
@test mean.(sysh.C) ≈ mean.(sys.C) atol=1e-5
@test std.(sysh.C) ≈ std.(sys.C) rtol=1e-5
@test mean.(sysh.D) ≈ mean.(sys.D) atol=1e-5
@test std.(sysh.D) ≈ std.(sys.D) rtol=1e-5


## TODO: add example from p. 173 Zhou

δm = δ()
δc = δ()
δk = δ()

m̄ = 2
c̄ = 1
k̄ = 10

m = m̄*(1 + 0.1*δm)
c = c̄*(1 + 0.2*δc)
k = k̄*(1 + 0.3*δk)

M = [
    0 1 0 0 0 0
    -k̄/m̄ -c̄/m̄ 1/m̄ -1/m̄ -1/m̄ -0.1/m̄
    0.3k̄ 0 0 0 0 0
    0 0.2c̄ 0 0 0 0
    -k̄ -c̄ 1 -1 -1 -0.1
]

A = [
    0 1
    -k/m -c/m
]
B = [
    0
    1/m
]

sys = ss(A, B, I(2), 0)
using Optim, ComponentArrays
# Optim.similar_axis(x::ComponentArray, n) = x .* zeros(length(x), n)
unsafe_comparisons()
l, res = find_lft(sys, [δm, δc, δk])
# l, res = find_lft(sys, 3)
Mh = Matrix(l)
norm(M - Mh) / norm(M)
sysh = ss(l,sys)

@test mean.(sysh.A) ≈ mean.(sys.A) rtol=1e-2
@test std.(sysh.A) ≈ std.(sys.A) rtol=1e-2
@test mean.(sysh.B) ≈ mean.(sys.B) rtol=1e-2
@test std.(sysh.B) ≈ std.(sys.B) rtol=1e-2
