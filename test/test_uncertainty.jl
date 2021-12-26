using RobustAndOptimalControl, ControlSystems, MonteCarloMeasurements

d = δr()
@test d.val == 0
@test d.radius == 1

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
    @test size(p.M) == (i+j+k, i+j+k)
end


##
@test uss(ss(I(2))).M.D == I
@test sum(δss(2, 2).M.D) == 4
@test sum(δss(4, 4).M.D) == 8

##
H = δss(2, 3)
W = ss(tf([1, .1],[.1, 1])*I(2))
WH = W*H

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

temp = (ss(I(1)) + W*d)
@test temp.nu == temp.ny == 1
@test temp.nz == temp.nw == 1

@test length(d.Δ) == 1


Pd = Pn*(ss(I(1)) + W*d)


usyss = sminreal(system_mapping(Pd))
@test !iszero(usyss.B)
@test !iszero(usyss.C)
@test (usyss.C*usyss.B)[] == 60
@test usyss == Pn


Pp = rand(Pd, 50)
@test Pp.nx == 1+1+2 # == nom, W, 2 sample unc

if isinteractive()
    bodeplot(Pp, w)
    bodeplot!(Pn, w)
end

Pn = ssrand(3,4,5)
@test δss(4,4, bound=0.2).M.D == [0I(4) sqrt(0.2)*I(4); sqrt(0.2)*I(4) 0I(4)]
mu = ss(I(4)) + δss(4,4, bound=0.2)


@test mu.M.D == [0I(4) sqrt(0.2)*I(4); sqrt(0.2)*I(4) I(4)]
@test size(Pn) == size(convert(UncertainSS, Pn))
@test convert(UncertainSS, Pn).D22 == Pn.D



P = Pn*mu


Pn2 = system_mapping(P)
@test Pn2 == Pn
@test size(P.M) == (7,8)

