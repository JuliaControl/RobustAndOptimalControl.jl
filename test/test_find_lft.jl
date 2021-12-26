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

@test pmean.(sysh.A) ≈ pmean.(sys.A) atol=1e-5
@test pstd.(sysh.A) ≈ pstd.(sys.A) rtol=1e-5
@test pmean.(sysh.B) ≈ pmean.(sys.B) atol=1e-5
@test pstd.(sysh.B) ≈ pstd.(sys.B) rtol=1e-5
@test pmean.(sysh.C) ≈ pmean.(sys.C) atol=1e-5
@test pstd.(sysh.C) ≈ pstd.(sys.C) rtol=1e-5
@test pmean.(sysh.D) ≈ pmean.(sys.D) atol=1e-5
@test pstd.(sysh.D) ≈ pstd.(sys.D) rtol=1e-5


## TODO: add example from p. 173 Zhou

δm = δr(32)
δC = δr(32)
δk = δr(32)

m̄ = 2
c̄ = 1
k̄ = 10

m = m̄*(1 + 0.1*δm)
c = c̄*(1 + 0.2*δC)
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
l, res = find_lft(sys, [δm, δC, δk])
# l, res = find_lft(sys, 3)
Mh = Matrix(l)
norm(M - Mh) / norm(M)
sysh = ss(l,sys)

@test pmean.(sysh.A) ≈ pmean.(sys.A) rtol=1e-2
@test pstd.(sysh.A) ≈ pstd.(sys.A) rtol=1e-2
@test pmean.(sysh.B) ≈ pmean.(sys.B) rtol=1e-2
@test pstd.(sysh.B) ≈ pstd.(sys.B) rtol=1e-2
