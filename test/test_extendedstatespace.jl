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

