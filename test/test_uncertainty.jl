using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements, Test, LinearAlgebra

d = δr()
@test d.val == 0
@test d.radius == 1
@test size(d) == (1,)
@test size(d, 1) == 1
@test length(d) == 1
@test eltype(d) == Float64
@test promote(1, δr(big(1.0))) isa Tuple{δ{BigFloat, BigFloat}, δ{BigFloat, BigFloat}}
@test convert(δ{Float64, Float64}, 1).val == 1

d2 = δr(1.0)
@test (-d2).val == -1

@test δr(10) isa Particles{Float64, 10}
@test δc(10) isa Complex{Particles{Float64, 10}}



d2 = d+d
@test d2.val == 0
@test d2.radius == 2
@test d2.name === d.name

d0 = d - d
@test d0.val == 0
@test d0.radius == 0
@test d2.name === d.name

d_ = δr()
d0_ = d - d_
@test d0_ isa UncertainSS


sys = uss(d)
@test sys.nx == 0
@test sys.D == [0 1; d.radius d.val] # lft representation of uncertain scalar


D = RobustAndOptimalControl.Δ(2, δc)
@test rand(D, 100) isa Matrix{Complex{Particles{Float64, 100}}}



# Example from Mu tools user guide
d = 2.4 + 0.4δr()
sys = uss(d)
@test sys.D == [0 1; 0.4 2.4]
isys = inv(sys)
@test isys.D ≈ [-1/6 1/2.4; -1/6 1/2.4] # Mu user guide seems wrong here


P = uss([δc(), δc()])
@test iszero(P.D11)
@test P.D21 == I
@test P.D12 == I
@test iszero(P.D22)
# W = tf(1.0, [1, 1])
# W*P



δel = δr(2, 1, :δel)
η = δr(6, 1, :η)
ρ = δr(-1, 1, :ρ)

s1 = uss(δel)
s2 = uss(η)

@test (1/η).D ≈ [
    -0.1667    0.1667
    -0.1667    0.1667
] atol=1e-3


@test (δel*η).D ≈ [
    0     1     6
    0     0     1
    1     2    12
] atol = 1e-3

@test (δel*η + ρ).D == [
    0     1     0     6
    0     0     0     1
    0     0     0     1
    1     2     1    11
] 

@test (δel/η).D ≈ [
    0   -0.1667    0.1667
    0   -0.1667    0.1667
1.0000   -0.3333    0.3333
] atol=1e-3



@test svdvals([η+3+δel δel/η].D)[1] ≈ 11.1896 rtol=1e-2

A = [[3+δel+η δel/η]; [7+ρ ρ+δel*η]]
@test svdvals(A.D)[1] ≈ 15.3545 rtol=0.01

# lft(A.D, A.Δ) not yet supported
# A.Δ
# A2 = lft(ss(A.sys), ss(Diagonal(A.Δ)), :u)
# @test A2 ≈ A



δel_mcm = rand(δel, 100)
η_mcm = rand(η, 100)
ρ_mcm = rand(ρ, 100)

A_mcm = [[3+δel_mcm+η_mcm δel_mcm/η_mcm]; [7+ρ_mcm ρ_mcm+δel_mcm*η_mcm]]


As = rand(A, 100)
@test all(As.D .≈ A_mcm)



##
for i = 1:2, j = 1:2, k = 1:2
    d1 = δss(i,j)
    d2 = δss(j,k)
    p = d1*d2
    @test p.ny == d1.ny
    @test p.nu == d2.nu
    @test size(p.M) == (j+k, i+j)
end


##
@test uss(I(2)).M.D == I
@test sum(δss(2, 2).D) == 4
@test sum(δss(4, 4).D) == 8

##
H = δss(2, 3)
W = tf([1, .1],[.1, 1]) .* I(2)
WH = W*H
@test WH.D == [zeros(3,2) I(3); 10I(2) zeros(2,3)]

##

w = 2π .* exp10.(LinRange(-2, 2, 500))
W = makeweight(0.40,15,3)
Pn = tf(1, [1/60, 1]) |> ss
d = δss(1,1)

Gu = rand(d, 100)
@test Gu.nx == d.Δ[1].nx
@test Gu.ny == d.Δ[1].ny
@test Gu.nu == d.Δ[1].nu


@test d.ny == 1
@test d.nu == 1
@test d.nz == 1
@test d.nw == 1

@test d.zinds == 1:1
@test d.winds == 1:1
@test d.uinds == 2:2
@test d.yinds == 2:2


temp = (W*d)
@test temp.nu == temp.ny == 1
@test temp.nz == temp.nw == 1

temp = I(1) + W*d
@test temp.nu == temp.ny == 1
@test temp.nz == temp.nw == 1

@test length(d.Δ) == 1


Pd = Pn*(I(1) + W*d)


usyss = sminreal(system_mapping(Pd))
@test !iszero(usyss.B)
@test !iszero(usyss.C)
@test (usyss.C*usyss.B)[] == 60
@test usyss == Pn


Pp = rand(Pd, 200)
@test Pp.nx == 1+1+2 # == nom, W, 2 sample unc

Gcl = lft(Pd, ss(-1))
structured_singular_value(Gcl)
unsafe_comparisons(true)
mvnyquistplot(Pp, w, points=true, Ms_circles=[1.2, 1.5])

if isinteractive()
    bodeplot(Pp, w, ylims=(1e-1, 1e1))
    bodeplot!(Pn, w)
end
unsafe_comparisons(false)

Pn = ssrand(3,4,5)
@test δss(4,4, bound=0.2).D == [0I(4) sqrt(0.2)*I(4); sqrt(0.2)*I(4) 0I(4)]
mu = I(4) + δss(4,4, bound=0.2)


@test mu.D == [0I(4) sqrt(0.2)*I(4); sqrt(0.2)*I(4) I(4)]
@test size(Pn) == size(convert(UncertainSS, Pn))
@test convert(UncertainSS, Pn).D22 == Pn.D

P = Pn*mu
Pn2 = system_mapping(P)
@test Pn2 == Pn
@test size(P) == (7,8)
@test size(P.M) == (4,4)


## this time mimo real
delta = uss([δr(), δr()])
a = 1
P = ss([0 a; -a -1], I(2), [1 a; 0 1], 0)* (I(2) + delta)
K = ss(I(2))

G = lft(P, -K)
hn = norm(G, Inf)

w = exp10.(LinRange(-5, 1, 100))
M = freqresp(G.M, w)
# mu = mussv_DG(M)
# maximum(mu)
# # maximum(structured_singular_value(M))
# @test 1/maximum(mu) ≈ √(2) atol=0.01



## Test block structure
d = δr()
@test RobustAndOptimalControl.makeblock([d]) == [-1, 0]
@test RobustAndOptimalControl.makeblock([d,d]) == [-2, 0]

d = δc()
@test RobustAndOptimalControl.makeblock([d]) == [1,0]
@test RobustAndOptimalControl.makeblock([d,d]) == [2, 0]
@test RobustAndOptimalControl.makeblock(δss(2,3).Δ) == [2, 3]



Δ = [δr(0,1,:qkarne), δr(0,1,:chuck_norris)]
@test RobustAndOptimalControl.block_structure(Δ)[1] == [[-1, 0], [-1, 0]]

Δ = [δr(0,1,:qkarne), δc(0,1,:chuck_norris)]
@test RobustAndOptimalControl.block_structure(Δ)[1] == [[1, 0], [-1, 0]] # The two elements should have been sorted and 1,-1 are therefore switched


Δ = [δr(0,1,:a), δss(2,3).Δ[], δr(0,1,:a), δc(0,1,:chuck_norris), δr(0,1,:a)]

@test RobustAndOptimalControl.block_structure(Δ)[1] == [
    [2, 3], # The gensym name will come first
    [-3, 0],
    [1, 0],
]

perm = RobustAndOptimalControl.block_structure(Δ)[2]
@test perm == [2,1,3,5,4]

P = partition(ssrand(8, 7, 2), 6, 7)
display(P)
P = UncertainSS(P, Δ)
display(P)

blocks, M = RobustAndOptimalControl.blocksort(P)
@test blocks == [
    [2, 3], # The gensym name will come first
    [-3, 0],
    [1, 0],
]

@test M[7,6] == P.M[6, 5] # test that the permutation was correctly applied, only test for the complex permutation
@test M[1:3,1:2] == P.M[2:4, 2:3]# test that the permutation was correctly applied, only test for the matrix block


## Test structured_singular_value with different kinds of blocks
w = exp10.(LinRange(-2, 2, 500))
delta = δss(1,1)
P = (tf(1,[1, .2, 1])) * (1+0.2*delta)
s = tf("s")
K = ss(1 + 2/s + 0.9s/(0.1s+1))
Gcl = lft(P, -K)
muplot(Gcl, w) # to test that the plot not errors

# These should all be equivalent
# Full block complex uncertainty gives μ = hinfnorm
@test norm(Gcl, Inf) ≈ 1.1861236982687908 rtol=1e-3
mu = structured_singular_value(Gcl, w)
@test 1/norm(mu, Inf) ≈ 1/1.1861236982687908 atol=0.005
@test RobustAndOptimalControl.robstab(Gcl) ≈ 0.84 atol=0.005
@test 1/structured_singular_value(Gcl) ≈ 0.84 atol=0.005


## same as above but with scalar instead of 1×1 system
w = exp10.(LinRange(-2, 2, 500))
delta = δc()
P = (tf(1,[1, .2, 1])) * (1+0.2*delta)
s = tf("s")
K = ss(1 + 2/s + 0.9s/(0.1s+1))
Gcl = lft(P, -K)


# These should all be equivalent
# Full block complex uncertainty gives μ = hinfnorm
@test norm(Gcl, Inf) ≈ 1.1861236982687908 rtol=1e-3
mu = structured_singular_value(Gcl, w)
@test 1/norm(mu, Inf) ≈ 1/1.1861236982687908 atol=0.005
@test RobustAndOptimalControl.robstab(Gcl) ≈ 0.84 atol=0.005
@test 1/structured_singular_value(Gcl) ≈ 0.84 atol=0.005


## this time mimo complex
delta = uss([δc(), δc()])
P = ss([0 1; 0 0], I(2), [1 0], 0) * (I(2) + delta)
# diagonal input uncertainty

K = ss([1;1])
G = lft(P, -K)
hn = norm(G, Inf)
@test robstab(G) ≈ 0.4824 atol=0.001
dm = diskmargin(system_mapping(P), K, 0, w)
dmm = argmin(dm->dm.margin, dm.simultaneous_input)
@test dmm.margin ≈ 0.5335 atol = 0.001


## Uncertain delays
using RobustAndOptimalControl, ControlSystemsBase, MonteCarloMeasurements
unsafe_comparisons(true)
L = Particles(collect((1:4) ./ 1000)) # Uncertain time delay, an integer number of milliseconds between 1ms and 4ms
P = delay(L)*tf(1, [0.01, 1])
C = pid(2, 1)
w = exp10.(-1:0.01:4)

bodeplot(P, exp10.(-1:0.001:3))
nyquistplot(P*C, w[1:10:end], points=true, xlims=(-3.5, 2.5), ylims=(-5, 1.5), Ms_circles=[1.5, 2], alpha=1) # Note, the nyquistplot with uncertain coefficients requires manual selection of plot limits
# plot(step(feedback(P, C), 0:0.0001:0.05), lab="L = " .* string.(P.Tau[].particles'), title="Disturbance response") # This is not being run to avoid having to load ControlSystems



## 

P = ss(-1, 1, 1, 0)
Pd = P * delay(0.1 .. 0.3)
@test Pd.P.P ≈ (tf(P) * delay(0.1 .. 0.3)).P.P


## Convert a named particle system to a single system with multiple outputs
# Tests that output and state names are set in a reasonable way
N = 20
A = randn(2,2) .+ 0.1Particles(N)
B = randn(2,3) .+ 0.1Particles(N)
C = randn(2,2) .+ 0.1Particles(N)
D = 0
P = named_ss(ss(A,B,C,D), "P")

MOP = RobustAndOptimalControl.mo_sys_from_particles(P);

@test MOP.ny == P.ny*N
@test MOP.nx == P.nx*N
@test MOP.nu == P.nu

@test MOP.x[1] == :Px1_1
@test MOP.x[2] == :Px2_1
@test MOP.x[3] == :Px1_2
@test MOP.u == P.u

# res = step(MOP, 0:0.1:1, method=:tustin)

N = 2000
A = randn(2,2) .+ 0.1Particles(N)
B = randn(2,3) .+ 0.1Particles(N)
C = randn(2,2) .+ 0.1Particles(N)
D = 0
P = named_ss(ss(A,B,C,D), "P")

MOP = RobustAndOptimalControl.mo_sys_from_particles(P);

using SparseArrays
@test MOP.ny == P.ny*N
@test MOP.nx == P.nx*N
@test MOP.nu == P.nu

@test MOP.A isa SparseMatrixCSC


## Poles and zeros of tf

ω = 4..12
ζ = 0.054..0.084

P  = tf([2*ζ/ω, 1],  [1/ω^2, 2*ζ/ω, 1, 0, 0])

C = 1.0 * tf([1, 1], [0.1, 1])
G = feedback(C*P, 1)
poles(G)
tzeros(G)