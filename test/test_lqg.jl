using RobustAndOptimalControl, ControlSystemsBase, LinearAlgebra, Test
# NOTE: Glad, Ljung chap 9 contains several numerical examples that can be used as test cases


approxin(el,col;kwargs...) = any(colel -> isapprox(el, colel; kwargs...), col)
approxsetequal(s1,s2;kwargs...) = all(approxin(p,s1;kwargs...) for p in s2) && all(approxin(p,s2;kwargs...) for p in s1)


w = exp10.(range(-5, stop=5, length=1000))
s = tf("s")
P = [1/(s+1) 2/(s+3); 1/(s+1) 1/(s+1)]
sys = ss(P)
sysmin = minreal(sys)
A,B,C,D = sys.A,sys.B,sys.C,sys.D
Am,Bm,Cm,Dm = sysmin.A,sysmin.B,sysmin.C,sysmin.D

@test approxsetequal(eigvals(Am), [-3,-1,-1])

Q1 = 100ones(4)
Q2 = 1I(2)
R1 = 100I(4)
R2 = 1I(2)
G = LQGProblem(sys, Q1, Q2, R1, R2)
gangoffourplot(G) # Test that it at least does not error
@test approxsetequal(eigvals(observer_controller(G).A), [ -31.6209+0.0im, -1.40629+0.0im, -15.9993+0.911174im, -15.9993-0.911174im, ], rtol = 1e-3)

C = observer_controller(G)
Ce = extended_controller(G)
@test system_mapping(Ce) == -C

# Test some utility functions
C = observer_controller(G)
@test C.A == sys.A - sys.B*lqr(G) - kalman(G)*sys.C
@test C.C == lqr(G)
@test C.B == kalman(G)
@test all(iszero, C.D)
RobustAndOptimalControl.gangoffourplot(sys, C)
RobustAndOptimalControl.gangofsevenplot(sys, C, tf(1, [1, 1]) .* I(2))

Ce = extended_controller(G)
@test Ce.ny == C.ny
@test Ce.nu == C.nu
@test Ce.nw == size(Q1, 1)
@test system_mapping(Ce) == -C
Cp = performance_mapping(Ce)
@test size(Cp) == (0, 4)
@test all(iszero, Cp.B)
@test Cp.nu == Ce.nw

# testing with integral action
# qQ = 1
# qR = 1
# Q1 = 1000I(4)
# Q2 = 1I(2)
# R1 = 1I(6)
# R2 = 1I(2)
# N = I(6)
# Gi = LQGProblem(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true, ϵ=0.001, N=N)
# gangoffourplot(Gi) # Test that it at least does not error
# @test approxsetequal(eigvals(Gi.sysc.A), [-0.001, -0.001, -47.4832, -44.3442, -3.40255, -1.15355 ], rtol = 1e-2)

# @test approxsetequal(eigvals(closedloop(G).A), [-1.0, -14.1774, -2.21811, -14.3206, -1.60615, -22.526, -1.0, -14.1774], rtol=1e-3)
# @test approxsetequal(eigvals(G.T.A), [-22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)
# @test approxsetequal(eigvals(G.S.A), [-22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)
# @test approxsetequal(eigvals(G.CS.A), [-31.6209+0.0im, -1.40629+0.0im, -15.9993+0.911174im, -15.9993-0.911174im, -22.526+0.0im, -2.21811+0.0im, -1.60615+0.0im, -14.3206+0.0im, -14.1774+0.0im, -1.0+0.0im, -1.0+0.0im, -14.1774+0.0im], rtol=1e-3)
# @test approxsetequal(eigvals(G.PS.A), [-1.0, -1.0, -3.0, -1.0, -22.526, -2.21811, -1.60615, -14.3206, -14.1774, -1.0, -1.0, -14.1774], rtol=1e-3)


# @test approxsetequal(eigvals(Gi.cl.A), [-1.0, -44.7425, -44.8455, -2.23294, -4.28574, -2.06662, -0.109432, -1.31779, -0.78293, -1.0, -0.001, -0.001], rtol=1e-3)
# @test approxsetequal(eigvals(Gi.T.A), [-44.7425, -44.8455, -4.28574, -0.109432, -2.23294, -2.06662, -1.31779, -0.78293, -1.0, -1.0], rtol=1e-3)
# @test approxsetequal(eigvals(Gi.S.A), [-44.7425, -44.8455, -4.28574, -0.109432, -2.23294, -2.06662, -1.31779, -0.78293, -1.0, -1.0], rtol=1e-3)

# @test eigvals(Gi.sysc.A) == eigvals(Gi.Fr.A)


## Test cases from Glad Ljung chap 9

# Aircraft control

A = [-0.292 8.13 -201 9.77 0 -12.5 17.1
    -0.152  -2.54  0.561  -0.0004  0  107  7.68
    0.0364  -0.0678  -0.481  0.0012  0  4.67  -7.98
    0 1 0.0401 0 0 0 0
    0 0 1 0 0 0 0
    0 0 0 0 0 -20 0
    0 0 0 0 0 0 -20
]

B = [
    0 -2.15
    -31.7 0.0274
    0 1.48
    0 0
    0 0
    20 0
    0 20
]
# Example 9.1
M = [0 0 0 1 0 0 0
     0 0 0 0 1 0 0]
C = M

N = I(7)
nx,ny,nu = size(A,1), size(C, 1), size(B,2)


Q1 = I(ny)
Q2 = I(nu)
R1 = I(nx)
R2 = I(ny)


sys = ss(A,N,B,M,C)
G = LQGProblem(sys, Q1, Q2, R1, R2)
@test lqr(G) ≈ [-0.0022 0.17 0.12 0.98 0.31 0.76 0.018
            -0.0038 0.028  -0.25 0.16  -0.95 0.072 0.10] rtol=0.01


Q2 = 0.1*I(nu)
G = LQGProblem(sys, Q1, Q2, R1, R2)
@test lqr(G) ≈ [-0.0036 0.39 0.21 3.11 0.63 1.54 0.046
            -0.0073 0.085 -0.78 0.57 -3.10 0.20 0.30] rtol=0.01

@test dcgain(closedloop(G)*static_gain_compensation(G)) ≈ I



# Room temperature control
K1 = 1
K2 = 1
K3 = 1/8
K4 = 0.2
A = [
    -K3-K4 K4 0 0 0
    K1 -K1-K2 0 K2 0
    0 0 -0.01 -0.068 0
    0 0 1 0 0
    0 0 0 0 -0.0001
]
# @warn "I added the last -0.0001"

B2 = [
    K3
    0
    0
    0
    0
]

C1 = [0 1 0 0 0]
C2 = [
    0 1 0 0 1
    0 1 0 0 0
    0 0 0 1 0
]
B1 = [
    0 0
    0 0
    1 0
    0 0
    0 1
]
sys = ss(A,B1,B2,C1,C2)


Q1 = 1.0I(1)
Q2 = 0.01I(1)
R1 = diagm([10.0, 1])
# R1 = diagm([0, 0, 10.0, 1, 0])
R2 = diagm([0.0001, 0.1, 0.001])
G = LQGProblem(sys, Q1, Q2, R1, R2, qQ = 0.000001)

# NOTE
# the kalman filter qutions will not work out until the noise model is incorporated

@test kalman(G) ≈ [-0.0039  0.0003  0.0035
            0.0463  0.0027  0.9676
            8.5224  0.0148  99.7543
            1.1426  0.0097  14.1198
            99.9634  -0.0026 -0.8533] rtol=1e-2
@test lqr(G) ≈ [2.8106 1.4071 6.8099 3.7874 0] rtol=1e-2

@test dcgain(RobustAndOptimalControl.closedloop(G)*static_gain_compensation(G))[] ≈ 1


## example 9.3

A = diagm([-1.0, -1, -2])
B = [
    1.0 0
    0.0 1
    0.0 1
]

C = [
    2.0 0 3
    1.0 1 0
]
M = C
sys = ss(A,I(3), B, M, C)


Q1 = I(2)
Q2 = 0.1I(2)
R1 = I(2)
R2 = 0.001I(2)

# test with integrator
# G = LQGProblem(sys, Q1, Q2, R1, R2, integrator=true, M = M, measurement=true)

# Lr = pinv(M * ((B * lqr(G)[:, 1:3] - A) \ B))

# @test Lr ≈ [
#     4.24 -3.43
#     0.12 4.82
# ] rtol = 0.01

# # NOTE: Glad ljung does not have identity in the end here, but the eigenvalues of A-BL are identical for this L and their L
# @test lqr(G) ≈ [
#     4.05 0.81 4.23 1 0 
#     5.05 1.01 5.96 0 1
# ] rtol = 0.01

# @test kalman(G) ≈ [
#     0 0
#     0 0
#     0 0
#     31.62 0
#     0 31.62
# ] rtol = 0.01

# @test dcgain(closedloop(G)*static_gain_compensation(G)) ≈ I

# Random system, compare with ML output


A = [  0.154388  1.42524    0.914886  -0.313071
1.36189   0.30188   -0.767355   0.313058
-0.485766  0.926861  -0.755906  -0.309824
2.10783   0.146366  -1.65582    0.115643]

B = [
  0.625203  -1.87059
 -0.399363   0.516783
 -0.150984   0.161573
 -0.432394  -0.162178
]

C = [
  0.99517   -0.26967    0.255771  -0.0356366
 -0.358348   0.114812  -0.417828   0.983806
  0.45019   -0.166284   1.47825    2.21589
]

K = [ 1.22237  -0.114471  -1.10722
1.00941   0.445568  -1.91614]

P = ss(A,B,C,0)
C = ss(K)

SiA = [
 -0.3534     1.707    -3.451    -5.871
   1.674    0.1881    0.1495     1.257
 -0.3683    0.8871   -0.5026  -0.08282
   2.433   0.09358    -2.655    -1.636
]

SiB = [
  0.6252   -1.871
 -0.3994   0.5168
  -0.151   0.1616
 -0.4324  -0.1622
]

SiC = [
   -0.759    0.1587     1.276      2.61
  0.01776  -0.09757     2.761     3.844
]

SiD = [
    1      0
    0      1]

Si = ss(SiA, SiB, SiC, SiD)


@test isapprox(input_sensitivity(P, C), Si, rtol=1e-3)


SoA = [
    -0.3534     1.707    -3.451    -5.871
      1.674    0.1881    0.1495     1.257
    -0.3683    0.8871   -0.5026  -0.08282
      2.433   0.09358    -2.655    -1.636
]

SoB = [
    1.124     0.905    -2.892
 -0.03348    -0.276     0.548
  0.02146  -0.08928    0.1424
   0.6922   0.02276   -0.7895
]

SoC = [
   0.9952   -0.2697    0.2558  -0.03564
  -0.3583    0.1148   -0.4178    0.9838
   0.4502   -0.1663     1.478     2.216 
]

SoD = I(3)

So = ss(SoA, -SoB, -SoC, SoD)


@test isapprox(output_sensitivity(P, C), So, rtol=1e-3)



ATo  = [
      -0.3534     1.707    -3.451    -5.871
        1.674    0.1881    0.1495     1.257
      -0.3683    0.8871   -0.5026  -0.08282
        2.433   0.09358    -2.655    -1.636
]
 
BTo = [
    1.124     0.905    -2.892
    -0.03348    -0.276     0.548
    0.02146  -0.08928    0.1424
    0.6922   0.02276   -0.7895
    ]

CTo = [
    -0.9952   0.2697  -0.2558  0.03564
    0.3583  -0.1148   0.4178  -0.9838
    -0.4502   0.1663   -1.478   -2.216
]
To = ss(ATo, -BTo, -CTo, 0)

@test isapprox(output_comp_sensitivity(P, C), To, rtol=1e-3)

ATi = [
      -0.3534     1.707    -3.451    -5.871
        1.674    0.1881    0.1495     1.257
      -0.3683    0.8871   -0.5026  -0.08282
        2.433   0.09358    -2.655    -1.636
]
 
BTi = [
    0.6252   -1.871
    -0.3994   0.5168
    -0.151   0.1616
    -0.4324  -0.1622
    ]

CTi = [
0.759   -0.1587    -1.276     -2.61
    -0.01776   0.09757    -2.761    -3.844
]
 
Ti = ss(ATi, BTi, CTi, 0)
 
@test isapprox(input_comp_sensitivity(P, C), Ti, rtol=1e-3)


Gfbc = RobustAndOptimalControl.feedback_control(P, C)
@test linfnorm((Gfbc - [output_comp_sensitivity(P,C); G_CS(P,C)]))[1] < 1e-10 # numerical problems calculating a reduced model for the difference. They should be equal everywhere

## dare3  
using ControlSystemsBase, RobustAndOptimalControl, LinearAlgebra, Plots, Test
nx = 5
nu = 2
ny = 3
P = ssrand(ny,nu,nx, Ts = 1, proper=true)
Pd = add_input_differentiator(P)
Q1 = randn(nx, nx); Q1 = Q1*Q1' + 1e-3*I
Q2 = randn(nu, nu); Q2 = Q2*Q2' + 1e-3*I
Q3 = 10randn(nu, nu); Q3 = Q3*Q3' + 1e-3*I
x0 = randn(nx)
#
QN = dare3(P, Q1, Q2, Q3)
QNf = dare3(P, Q1, Q2, Q3, full=true)
L = lqr3(P, Q1, Q2, Q3)
Lf = lqr3(P, Q1, Q2, Q3, full=true)
res = lsim(P, (x,t)->-L*x, 5000; x0)
resf = lsim(Pd, (x,t)->-Lf*x, 5000; x0=[x0; zeros(nu)])

# actual_cost = dot(res.x, Q1, res.x) + dot(res.u, Q2, res.u) + 
# dot(res.u[:, 2:end] - res.u[:, 1:end-1], Q3, res.u[:, 2:end] - res.u[:, 1:end-1]) + dot(res.u[:, 1], Q3, res.u[:, 1])
# predicted_cost = dot(x0, QN, x0)
# @show (actual_cost - predicted_cost) / actual_cost
# NOTE: can't test this because the non-full L controller is not aware of the input state like an MPC controller with fixed u0 in the optimization cost would be.

Δu = resf.u[:, 2:end] - resf.u[:, 1:end-1]
actual_cost = dot(resf.x[1:end-nu, :], Q1, resf.x[1:end-nu, :]) +
              dot(resf.u, Q2, resf.u) +
              dot(Δu, Q3, Δu) +
              dot(resf.u[:, 1], Q3, resf.u[:, 1])

predicted_cost = dot(x0, QN, x0)
predicted_costf = dot([x0; zeros(nu)], QNf, [x0; zeros(nu)])
@test predicted_costf ≈ predicted_cost

# @show (actual_cost - predicted_cost) / actual_cost
@test actual_cost ≈ predicted_cost rtol=1e-10

## Double mass model with output penalty

# One input, to z outputs
P = DemoSystems.double_mass_model(outputs = 1:4)
# Pe = partition(P, y = 1:2, u = 1)
Pe = ExtendedStateSpace(P, C1 = P.C[3:4, :], C2 = P.C[1:1, :], B1 = P.C[[2,4], :]')

Q1 = 1.0I(Pe.nz)
Q2 = 0.1I(Pe.nu)
R1 = 1.0I(Pe.nw)
R2 = 0.1I(Pe.ny)

G = LQGProblem(Pe, Q1, Q2, R1, R2)
Cfb = observer_controller(G)
Ce = extended_controller(G)

@test Cfb.nu == Pe.ny
@test Cfb.ny == Pe.nu


Cff = RobustAndOptimalControl.ff_controller(G) # TODO: need to add argument that determins whether to use state references or references of z
@test Cff.nu == Pe.nz
@test Cff.ny == Pe.nu

Gcl2 = feedback(system_mapping(Pe), Cfb) * Cff
@test size(Gcl2) == (1, 2)
@test isstable(Gcl2)
@test dcgain(Gcl2) ≈ [1.0 0] atol=1e-6



## One input, one z output
P = DemoSystems.double_mass_model(outputs = 1:4)
# Pe = partition(P, y = 1:2, u = 1)
Pe = ExtendedStateSpace(P, C1 = P.C[3:3, :], C2 = P.C[1:1, :], B1 = P.C[[2,4], :]')

Q1 = 1.0I(Pe.nz)
Q2 = 0.1I(Pe.nu)
R1 = 1.0I(Pe.nw)
R2 = 0.1I(Pe.ny)

G = LQGProblem(Pe, Q1, Q2, R1, R2)
Cfb = observer_controller(G)
Ce = extended_controller(G)

@test Cfb.nu == Pe.ny
@test Cfb.ny == Pe.nu


Cff = RobustAndOptimalControl.ff_controller(G)
@test Cff.nu == Pe.nz
@test Cff.ny == Pe.nu

Gcl = feedback(system_mapping(Pe) * Cfb)
@test size(Gcl) == (1, 1)
@test isstable(Gcl)
@test dcgain(Gcl)[] ≈ 1.0

Gcl2 = feedback(system_mapping(Pe), Cfb) * Cff
@test size(Gcl2) == (1, 1)
@test isstable(Gcl2)
@test dcgain(Gcl2)[] ≈ 1.0

## Multibody cartpole tests

lsys = let
    lsysA = [0.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0; 4.474102070258828 0.0 0.0 0.0; 39.683507165892394 0.0 0.0 0.0]
    lsysB = [0.0; 0.0; 0.8415814526118872; 2.3385955754360115;;]
    lsysC = [0.0 1.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    lsysD = [0.0; 0.0; 0.0; 0.0;;]
    named_ss(ss(lsysA, lsysB, lsysC, lsysD), x=[:revolute₊phi, :prismatic₊s, :prismatic₊v, :revolute₊w], u=[Symbol("u(t)")], y=[Symbol("x(t)"), Symbol("phi(t)"), Symbol("v(t)"), Symbol("w(t)")])
end

C = lsys.C
Q1 = Diagonal([10, 10, 10, 1])
Q2 = Diagonal([0.1])

R1 = lsys.B*Diagonal([1])*lsys.B'
R2 = Diagonal([0.01, 0.01])
Pn = lsys[[:x, :phi], :]
P = ss(Pn)

lqg = LQGProblem(P, Q1, Q2, R1, R2)

# Below we test that the closed-loop DC gain from references to cart position is 1
Ce = extended_controller(lqg)


# Method 1: compute DC gain compensation using ff_controller, this is used as reference for the others. I feel uneasy about the (undocumented) comp_dc = false argument here, ideally a more generally applicable function ff_controller would be implemented that can handle both state and output references etc.
dc_gain_compensation = dcgain(RobustAndOptimalControl.ff_controller(lqg, comp_dc = false))[]

# Method 2, using the observer controller and static_gain_compensation
Cfb = observer_controller(lqg)
@test inv(dcgain((feedback(P, Cfb)*static_gain_compensation(lqg))[1,2]))[] ≈ dc_gain_compensation rtol=1e-8

# Method 3, using the observer controller and the ff_controller
cl = feedback(system_mapping(lqg), observer_controller(lqg))*RobustAndOptimalControl.ff_controller(lqg, comp_dc = true)
@test dcgain(cl)[1,2] ≈ 1 rtol=1e-8
@test isstable(cl)

# Method 4: Build compensation into R and compute the closed-loop DC gain, should be 1
R = named_ss(ss(dc_gain_compensation*I(4)), "R") # Reference filter
Ce = named_ss(ss(Ce); x = :xC, y = :u, u = [R.y; :y^lqg.ny])

Cry = RobustAndOptimalControl.connect([R, Ce]; u1 = R.y, y1 = R.y, w1 = [R.u; :y^lqg.ny], z1=[:u])

connections = [
    :u => :u
    [:x, :phi] .=> [:y1, :y2]
]
cl = RobustAndOptimalControl.connect([lsys, Cry], connections; w1 = R.u, z1 = [:x, :phi])
@test inv(dcgain(cl)[1,2]) ≈ 1 rtol=1e-8
@test isstable(cl)


# Method 5: close the loop manually with reference as input and position as output

R = named_ss(ss(I(4)), "R") # Reference filter, used for signal names only
Ce = named_ss(ss(extended_controller(lqg)); x = :xC, y = :u, u = [:Ry^4; :y^lqg.ny])
cl = feedback(lsys, Ce, z1 = [:x], z2=[], u2=:y^2, y1 = [:x, :phi], w2=[:Ry2], w1=[], pos_feedback=true)
@test inv(dcgain(cl)[]) ≈ dc_gain_compensation rtol=1e-8
@test isstable(cl)

cl = feedback(lsys, Ce, z1 = [:x, :phi], z2=[], u2=:y^2, y1 = [:x, :phi], w2=[:Ry2], w1=[], pos_feedback=true)
@test pinv(dcgain(cl)) ≈ [dc_gain_compensation 0] atol=1e-8
@test isstable(cl)

# Method 6: use the z argument to extended_controller to compute the closed-loop TF

Ce, cl = extended_controller(lqg, z=[1, 2])
@test pinv(dcgain(cl)[1,2]) ≈ dc_gain_compensation atol=1e-8
@test isstable(cl)

Ce, cl = extended_controller(lqg, z=[1])
@test pinv(dcgain(cl)[1,2]) ≈ dc_gain_compensation atol=1e-8
@test isstable(cl)

## Output references
Li = [1 0]
Cry2, cl = extended_controller(lqg; output_ref = true, Li, z=[1])
cl = minreal(cl)

@test dcgain(cl) ≈ [1 0] atol=1e-12
@test isstable(cl)

