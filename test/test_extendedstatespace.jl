using RobustAndOptimalControl
using ControlSystemsBase
using ControlSystemsBase: balance_statespace
using Test

G = ssrand(3,4,5)

G = partition(G, 1, 1)

@test G isa ExtendedStateSpace
@test G == G
@test G ≈ G
@test ControlSystemsBase.noutputs(G) == 3
@test ControlSystemsBase.ninputs(G) == 4
@test ndims(G) == 2
@test size(G) == (3,4)
@test eltype(G) == typeof(G)
@test ControlSystemsBase.numeric_type(G) == Float64
@test_throws ErrorException G[1]

@show partition(ssrand(2,2,1), 1, 1)

@test noise_mapping(G) == ss(G)[2:end, 1]

# Test promotion of regular system
G = ssrand(3,4,5)
Ge = ExtendedStateSpace(G)
@test performance_mapping(Ge) == G
@test system_mapping(Ge) == G


Ge = convert(ExtendedStateSpace{Continuous, Float64}, G)
@test system_mapping(Ge) == G
@test size(performance_mapping(Ge)) == (0,0)
@test ss(Ge) == G

Ge2, G2 = promote(Ge, G)
@test Ge2 == Ge
@test G2 == partition(G, 0, 0) # Promoted to ess with empty performance model.

## Connect and feedback_control

G = ssrand(3,4,2)
K = ssrand(4,3,2)

Gcl1 = feedback_control(G, K)
G = named_ss(G, :G)
K = named_ss(K, :K)
S = sumblock("Ku = r - Gy", n=3)

@test S[:, :r] == S[:, :r^3] # test indexing with a prefix

z1 = [G.y; K.y]
w1 = :r^3
connections = [K.y .=> G.u; G.y .=> G.y; K.u .=> K.u]
Gcl2 = connect([G, K, S], connections; z1, w1)

@test linfnorm(minreal(Gcl1 - Gcl2.sys))[1] < sqrt(eps())


##
G = ssrand(3,4,5)
Ge = partition(G, 1, 1)
@test ss(balance_statespace(Ge)[1]) == balance_statespace(G)[1]
@test ss(balreal(Ge)[1]) == balreal(G)[1]
@test ss(modal_form(Ge)[1]) == modal_form(G)[1]
@test ss(schur_form(Ge)[1]) == schur_form(G)[1]
@test ss(hess_form(Ge)[1]) == hess_form(G)[1]


## feedback
P = ssrand(2,3,2, proper=true)
C = ssrand(3,2,2)

@test feedback(P, C) == feedback(ExtendedStateSpace(P), C)
@test feedback(P, C) == feedback(P, convert(ExtendedStateSpace{Continuous, typeof(C)}, C))

@test feedback(P, C) == feedback(ExtendedStateSpace(P), convert(ExtendedStateSpace{Continuous, typeof(C)}, C))


P = ssrand(2,3,2, proper=false)
C = ssrand(3,2,2)

@test feedback(P, C) == feedback(ExtendedStateSpace(P), C)
@test feedback(P, C) == feedback(P, convert(ExtendedStateSpace{Continuous, typeof(C)}, C))

@test feedback(P, C) == feedback(ExtendedStateSpace(P), convert(ExtendedStateSpace{Continuous, typeof(C)}, C))