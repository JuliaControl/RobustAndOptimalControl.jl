using RobustAndOptimalControl, ControlSystems

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



dH = δss(2, 3)
W = ss(tf([1, .1],[.1, 1])*I(2))
WdH = W*dH

# s1 = W
# s2 = dH
# s1, s2 = promote(s1, s2)




# s1 = ssrand(4,4,4)
# s2 = ssrand(4,4,4)

# se1 = partition(s1,2,2)
# se2 = partition(s2,2,2)

# @test ss(se1*se2) ≈ (s1*s2) rtol=1e-10

Wunc = makeweight(0.40,15,3)
sysNom = tf(1, [1/60, 1]) |> ss
unc = δss(1,1)
@test unc.ny == 1
@test unc.nu == 1
@test unc.nz == 2
@test unc.nw == 2

@test unc.zinds == 1:2
@test unc.winds == 1:2
@test unc.uinds == 3:3
@test unc.yinds == 3:3


temp = (Wunc*unc)
@test temp.nu == temp.ny == 1
@test temp.nz == temp.nw == 2

temp = (ss(I(1)) + Wunc*unc)
@test temp.nu == temp.ny == 1
@test temp.nz == temp.nw == 2



usys = sysNom*(ss(I(1)) + Wunc*unc)


# @test usys.D11 == [0 1; 1 0]
usyss = system_mapping(usys)
@test !iszero(usyss.B)
@test !iszero(usyss.C)
@test (usyss.C*usyss.B)[] == 60