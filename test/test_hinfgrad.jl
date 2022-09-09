using RobustAndOptimalControl, ControlSystemsBase
using FiniteDiff, Zygote



##
# Define the proces parameters
k1, k2, kc, g = 3.33, 3.35, 0.5, 981
A1, A3, A2, A4 = 28, 28, 32, 32
a1, a3, a2, a4= 0.071, 0.071, 0.057, 0.057
h01, h02, h03, h04 = 12.4, 12.7, 1.8, 1.4
T1, T2 = (A1/a1)*sqrt(2*h01/g), (A2/a2)*sqrt(2*h02/g)
T3, T4 = (A3/a3)*sqrt(2*h03/g), (A4/a4)*sqrt(2*h04/g)
c1, c2 = (T1*k1*kc/A1), (T2*k2*kc/A2)
γ1, γ2 = 0.7, 0.6

# Define the process dynamics
A = [-1/T1     0 A3/(A1*T3)          0;
     0     -1/T2          0 A4/(A2*T4);
     0         0      -1/T3          0;
     0         0          0      -1/T4];
B = [γ1*k1/A1     0;
     0                γ2*k2/A2;
     0                (1-γ2)*k2/A3;
     (1-γ1)*k1/A4 0              ];
C = [kc 0 0 0;
     0 kc 0 0];
D = zeros(2,2)
G = ss(A,B,C,D);

# Sensitivity weight function
WS = makeweight(100, (0.1, 1), 1/10) * I(2)

# Output sensitivity weight function
# WUelement = 5*tf([1,1],[0.1,1])
WU = tf(0.01) .* I(2)

# Complementary sensitivity weight function
# WT = tf([10,0.1],[1,1]) .* I(2)
WT = makeweight(1/10, (0.1, 1), 10) * I(2)

# Form augmented P dynamics in state-space
P = hinfpartition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hinfassumptions(P)

# Synthesize the H-infinity optimal controller
C, γ = hinfsynthesize(P)

Pcl, S, CS, T = hinfsignals(P, G, C)
isinteractive() && specificationplot([S, CS, T], [WS[1,1], 0.01, WT[1,1]], γ)



function cost(v)
    K = vec2sys(v, 2, 2)
    T = lft(P, K)
    hinfnorm2(T)[1]
end

v = vec(C)

# v0 = 0.5*vec(C)
# cost(v0)

# K0 = vec2sys(v0, 2, 2)
# T0 = lft(P, K0)
# using Optim, Optim.LineSearches
# res = Optim.optimize(
#     cost,
#     v0,
#     ParticleSwarm(),
#     Optim.Options(
#         store_trace       = true,
#         show_trace        = true,
#         show_every        = 1,
#         iterations        = 1000,
#         allow_f_increases = false,
#         time_limit        = 1,
#         x_tol             = 0,
#         f_abstol          = 0,
#         g_tol             = 1e-8,
#         f_calls_limit     = 0,
#         g_calls_limit     = 0,
#     ),
# )

# vo = res.minimizer
# Ko = vec2sys(vo, 2, 2)
# To = lft(P, Ko)

# Pcl, S, CS, T = hinfsignals(P, G, C)
# specificationplot([S, CS, T], [WS[1,1], 0.01, WT[1,1]], 1, nsigma=1)

# Pcl, S, CS, T = hinfsignals(P, G, Ko)
# specificationplot!([S, CS, T], l=(:black, :dash), nsigma=1, ylims=(0.01, Inf))


##

function testfun(v)
    K = vec2sys(v, 2, 2)
    typeof(K)
    hinfnorm2(K, rtolinf=1e-8)[1]
end

testfun(v)

g1 = Zygote.gradient(testfun, v)[1]

g2 = FiniteDiff.finite_difference_gradient(testfun, v, absstep=1e-8)

senssys1 = vec2sys(g1, 2, 2)
senssys2 = vec2sys(g2, 2, 2)

@test norm(senssys1.A - senssys2.A) / norm(senssys1.A) < 0.001
@test norm(senssys1.B - senssys2.B) / norm(senssys1.B) < 0.001
@test norm(senssys1.C - senssys2.C) / norm(senssys1.C) < 0.001
@test norm(senssys1.D - senssys2.D) / norm(senssys1.D) < 0.001


## Broken due to Zygote not handling try/catch and @warn. The offending method is feedback
# g1 = Zygote.gradient(cost, v)[1]
# g2 = FiniteDiff.finite_difference_gradient(cost, v, absstep=1e-8)

# senssys1 = vec2sys(g1, 2, 2)
# senssys2 = vec2sys(g2, 2, 2)

# @test norm(senssys1.A - senssys2.A) / norm(senssys1.A) < 0.001
# @test norm(senssys1.B - senssys2.B) / norm(senssys1.B) < 0.001
# @test norm(senssys1.C - senssys2.C) / norm(senssys1.C) < 0.001
# @test norm(senssys1.D - senssys2.D) / norm(senssys1.D) < 0.001
