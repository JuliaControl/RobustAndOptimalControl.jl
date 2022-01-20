using RobustAndOptimalControl
using ControlSystems

G = ssrand(3,4,5)

G = partition(G, 1, 1)

@test G isa ExtendedStateSpace
@test G == G
@test G ≈ G
@test ControlSystems.noutputs(G) == 3
@test ControlSystems.ninputs(G) == 4
@test ndims(G) == 2
@test size(G) == (3,4)
@test eltype(G) == typeof(G)
@test ControlSystems.numeric_type(G) == Float64
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

