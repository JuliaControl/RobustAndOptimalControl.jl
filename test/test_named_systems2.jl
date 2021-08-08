using RobustAndOptimalControl, ControlSystems
using RobustAndOptimalControl: @check_unique, @check_all_unique

G1 = ss(1,1,1,1)
G2 = ss(1,1,1,1)
s1 = named_ss(G1, x = [:x], u = :u, y=:y) # test providing only symbol
s2 = named_ss(G2, x = [:z], u = [:u], y=[:y])

@test s1[:y, :u] == s1
@test s1[[:y], [:u]] == s1

G3 = ControlSystems.ssrand(1,2,3)
s3 = named_ss(G3, x = [:x1, :x2, :x3], u = [:u1, :u2], y=[:y])

@test s3[:y, :u1] == named_ss(G3[1,1], x = [:x1, :x2, :x3], u = [:u1], y=[:y])
@test s3[:y, :u2] == named_ss(G3[1,2], x = [:x1, :x2, :x3], u = [:u2], y=[:y])
@test s3[1, 2] == named_ss(G3[1,2], x = [:x1, :x2, :x3], u = [:u2], y=[:y]) # indexing with integer preserves names

@test s3[:y, [:u1, :u2]] == s3

for op in (+, -)
    @testset "$op" begin
        @info "Testing $op"
        
        G12 = op(G1, G2)
        s12 = op(s1, s2)

        @test s12.x == [:x, :z]
        @test s12.u == [:u]
        @test s12.y == [:y]
    
        @test s12.A == G12.A
        @test s12.B == G12.B
        @test s12.C == G12.C
        @test s12.D == G12.D
    end
end


s1 = named_ss(G1, x = [:x], u = [:u1], y=[:y1])
s2 = named_ss(G2, x = [:z], u = [:u2], y=[:y2])

@testset "*" begin
    G12 = *(G1, G2)
    s12 = *(s1, s2)

    @test s12.x == [:x, :z]
    @test s12.u == [:u2]
    @test s12.y == [:y1]

    @test s12.A == G12.A
    @test s12.B == G12.B
    @test s12.C == G12.C
    @test s12.D == G12.D

    @test_throws ArgumentError s1*s1
end

@testset "Concatenation" begin
    @info "Testing Concatenation"
    
    # same output name
    s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
    s2 = named_ss(G2, x = :z, u = :u2, y=:y1)
    s3 = named_ss(G2, x = :w, u = :u3, y=:y1)
    s1s2 = [s1 s2]
    @test s1s2.sys == [G1 G2]

    @test s1s2.x == [:x, :z]
    @test s1s2.u == [:u1, :u2]
    @test s1s2.y == [:y1]

    @test_throws ArgumentError [s1 s2 s1 s2]
    @test_nowarn [s1 s2 s3]
    
    # same input name
    s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
    s2 = named_ss(G2, x = :z, u = :u1, y=:y2)
    s3 = named_ss(G2, x = :w, u = :u1, y=:y3)
    s1s2 = [s1; s2]
    @test s1s2.sys == [G1; G2]

    @test s1s2.x == [:x, :z]
    @test s1s2.u == [:u1]
    @test s1s2.y == [:y1, :y2]

    @test_throws ArgumentError [s1; s2; s1; s2]
    @test_nowarn [s1; s2; s3]
end

G1 = ss(1,1,1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
s2 = named_ss(G2, x = :z, u = :u2, y=:y2)
@test_nowarn @check_all_unique s1 s2

fb = feedback(s1, s2, r = :r)
# @test fb.x == [:x, :z]
@test fb.u == [:r]
@test fb.y == [:y1]


s1 = named_ss(G1, x = [:x], u = [:u1], y=[:y1])
s2 = named_ss(G2, x = [:z], u = [:u1], y=[:y2])
@test_throws ArgumentError @check_all_unique s1 s2

## Promotion and conversion
@testset "Promotion and conversion" begin
    @info "Testing Promotion and conversion"
    # Named and not named
    G1 = named_ss(ssrand(1,1,1))
    G2 = HeteroStateSpace(ssrand(1,1,2))
    @test promote_type(typeof.((G1, G2))...) == NamedStateSpace{Continuous, HeteroStateSpace{Continuous, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}

    G1p, G2p = promote(G1, G2)
    @test G1p isa NamedStateSpace{Continuous, HeteroStateSpace{Continuous, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}
    @test G2p isa NamedStateSpace{Continuous, HeteroStateSpace{Continuous, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}

    # Named and named
    @test promote_type(typeof.((G1, G2p))...) == NamedStateSpace{Continuous, HeteroStateSpace{Continuous, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}
    G1p, G2p2 = promote(G1, G2p)
    @test G1p isa NamedStateSpace{Continuous, HeteroStateSpace{Continuous, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}
    @test G2p == G2p
end