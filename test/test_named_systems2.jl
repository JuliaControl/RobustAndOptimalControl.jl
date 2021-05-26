using RobustAndOptimalControl, ControlSystems
using RobustAndOptimalControl: @check_unique, @check_all_unique

G1 = ss(1,1,1,1)
G2 = ss(1,1,1,1)
s1 = named_ss(G1, x_names = [:x], u_names = [:u], y_names=[:y])
s2 = named_ss(G2, x_names = [:z], u_names = [:u], y_names=[:y])

@test s1[:y, :u] == s1
@test s1[[:y], [:u]] == s1

G3 = ControlSystems.ssrand(1,2,3)
s3 = named_ss(G3, x_names = [:x1, :x2, :x3], u_names = [:u1, :u2], y_names=[:y])

@test s3[:y, :u1] == named_ss(G3[1,1], x_names = [:x1, :x2, :x3], u_names = [:u1], y_names=[:y])
@test s3[:y, :u2] == named_ss(G3[1,2], x_names = [:x1, :x2, :x3], u_names = [:u2], y_names=[:y])

@test s3[:y, [:u1, :u2]] == s3

for op in (+, -)
    @testset "$op" begin
        @info "Testing $op"
        
        G12 = op(G1, G2)
        s12 = op(s1, s2)

        @test s12.x_names == [:x, :z]
        @test s12.u_names == [:u]
        @test s12.y_names == [:y]
    
        @test s12.A == G12.A
        @test s12.B == G12.B
        @test s12.C == G12.C
        @test s12.D == G12.D
    end
end


s1 = named_ss(G1, x_names = [:x], u_names = [:u1], y_names=[:y1])
s2 = named_ss(G2, x_names = [:z], u_names = [:u2], y_names=[:y2])

@testset "*" begin
    G12 = *(G1, G2)
    s12 = *(s1, s2)

    @test s12.x_names == [:x, :z]
    @test s12.u_names == [:u2]
    @test s12.y_names == [:y1]

    @test s12.A == G12.A
    @test s12.B == G12.B
    @test s12.C == G12.C
    @test s12.D == G12.D

    @test_throws ArgumentError s1*s1
end

@testset "Concatenation" begin
    @info "Testing Concatenation"
    
    # same output name
    s1 = named_ss(G1, x_names = [:x], u_names = [:u1], y_names=[:y1])
    s2 = named_ss(G2, x_names = [:z], u_names = [:u2], y_names=[:y1])
    s3 = named_ss(G2, x_names = [:w], u_names = [:u3], y_names=[:y1])
    s1s2 = [s1 s2]
    @test s1s2.sys == [G1 G2]

    @test s1s2.x_names == [:x, :z]
    @test s1s2.u_names == [:u1, :u2]
    @test s1s2.y_names == [:y1]

    @test_throws ArgumentError [s1 s2 s1 s2]
    @test_nowarn [s1 s2 s3]
    
    # same input name
    s1 = named_ss(G1, x_names = [:x], u_names = [:u1], y_names=[:y1])
    s2 = named_ss(G2, x_names = [:z], u_names = [:u1], y_names=[:y2])
    s3 = named_ss(G2, x_names = [:w], u_names = [:u1], y_names=[:y3])
    s1s2 = [s1; s2]
    @test s1s2.sys == [G1; G2]

    @test s1s2.x_names == [:x, :z]
    @test s1s2.u_names == [:u1]
    @test s1s2.y_names == [:y1, :y2]

    @test_throws ArgumentError [s1; s2; s1; s2]
    @test_nowarn [s1; s2; s3]
end

G1 = ss(1,1,1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x_names = [:x], u_names = [:u1], y_names=[:y1])
s2 = named_ss(G2, x_names = [:z], u_names = [:u2], y_names=[:y2])
@test_nowarn @check_all_unique s1 s2

fb = feedback(s1, s2, r_names = :r)
@test fb.x_names == [:x, :z]
@test fb.u_names == [:r]
@test fb.y_names == [:y1]


s1 = named_ss(G1, x_names = [:x], u_names = [:u1], y_names=[:y1])
s2 = named_ss(G2, x_names = [:z], u_names = [:u1], y_names=[:y2])
@test_throws ArgumentError @check_all_unique s1 s2