using RobustAndOptimalControl, ComponentArrays, ControlSystems

G1 = ss(1,1,1,1)
G2 = ss(1,1,1,1)
s1 = named_ss(G1, state_names = [:x], input_names = [:u], output_names=[:y])
s2 = named_ss(G2, state_names = [:z], input_names = [:u], output_names=[:y])

@test s1[:y, :u] == s1
@test s1[[:y], [:u]] == s1

G3 = ControlSystems.ssrand(1,2,3)
s3 = named_ss(G3, state_names = [:x1, :x2, :x3], input_names = [:u1, :u2], output_names=[:y])

@test_broken s3[:y, :u1] == named_ss(G3[1,1], state_names = [:x1, :x2, :x3], input_names = [:u1], output_names=[:y])
@test_broken s3[:y, :u2] == named_ss(G3[1,2], state_names = [:x1, :x2, :x3], input_names = [:u2], output_names=[:y])

@test s3[:y, [:u1, :u2]] == s3

for op in (+, -)
    @testset "$op" begin
        @info "Testing $op"
        
        G12 = op(G1, G2)
        s12 = op(s1, s2)

        @test s12.A isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(x = 1, z = 2)}}}
        @test s12.B isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(u = 1,)}}}
        @test s12.C isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y = 1,)}, Axis{(x = 1, z = 2)}}}
        @test s12.D isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y = 1,)}, Axis{(u = 1,)}}}
        @test s12.A == G12.A
        @test s12.B == G12.B
        @test s12.C == G12.C
        @test s12.D == G12.D
    end
end


s1 = named_ss(G1, state_names = [:x], input_names = [:u1], output_names=[:y1])
s2 = named_ss(G2, state_names = [:z], input_names = [:u2], output_names=[:y2])

@testset "*" begin
    G12 = *(G1, G2)
    s12 = *(s1, s2)

    @test s12.A isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(x = 1, z = 2)}}}
    @test s12.B isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(u2 = 1,)}}}
    @test s12.C isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y1 = 1,)}, Axis{(x = 1, z = 2)}}}
    @test s12.D isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y1 = 1,)}, Axis{(u2 = 1,)}}}
    @test s12.A == G12.A
    @test s12.B == G12.B
    @test s12.C == G12.C
    @test s12.D == G12.D
end

@testset "Concatenation" begin
    @info "Testing Concatenation"
    
    # same output name
    s1 = named_ss(G1, state_names = [:x], input_names = [:u1], output_names=[:y1])
    s2 = named_ss(G2, state_names = [:z], input_names = [:u2], output_names=[:y1])
    s3 = named_ss(G2, state_names = [:w], input_names = [:u3], output_names=[:y1])
    s1s2 = [s1 s2]
    @test s1s2 == [G1 G2]
    @test s1s2.A isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(x = 1, z = 2)}}}
    @test s1s2.B isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(u1 = 1, u2 = 2)}}}
    @test s1s2.C isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y1 = 1,)}, Axis{(x = 1, z = 2)}}}
    @test s1s2.D isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y1 = 1,)}, Axis{(u1 = 1, u2 = 2)}}}

    @test_throws ArgumentError [s1 s2 s1 s2]
    @test_nowarn [s1 s2 s3]
    
    # same input name
    s1 = named_ss(G1, state_names = [:x], input_names = [:u1], output_names=[:y1])
    s2 = named_ss(G2, state_names = [:z], input_names = [:u1], output_names=[:y2])
    s3 = named_ss(G2, state_names = [:w], input_names = [:u1], output_names=[:y3])
    s1s2 = [s1; s2]
    @test s1s2 == [G1; G2]
    @test s1s2.A isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(x = 1, z = 2)}}}
    @test s1s2.B isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(x = 1, z = 2)}, Axis{(u1 = 1,)}}}
    @test s1s2.C isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y1 = 1, y2 = 2)}, Axis{(x = 1, z = 2)}}}
    @test s1s2.D isa ComponentArray{Int64,2,Matrix{Int64},Tuple{Axis{(y1 = 1, y2 = 2)}, Axis{(u1 = 1,)}}}

    @test_throws ArgumentError [s1; s2; s1; s2]
    @test_nowarn [s1; s2; s3]
end