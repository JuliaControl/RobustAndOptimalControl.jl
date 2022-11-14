using RobustAndOptimalControl, ControlSystemsBase
using RobustAndOptimalControl: @check_unique, @check_all_unique

@test :x^3 == expand_symbol(:x, 3) == [:x1, :x2, :x3]

G1 = ss(1.0,1,1,1)
G2 = ss(1.0,1,1,1)
s1 = named_ss(G1, x = [:x], u = :u, y=:y) # test providing only symbol
s2 = named_ss(G2, x = [:z], u = [:u], y=[:y])
@show s1
@test :x ∈ propertynames(s1)
@test :A ∈ propertynames(s1)

s_autoname = named_ss(G1, :G)
@test s_autoname.x == [:Gx]
@test s_autoname.y == [:Gy]
@test s_autoname.u == [:Gu]

@test ControlSystemsBase.state_names(s_autoname) == string.(s_autoname.x)
@test ControlSystemsBase.input_names(s_autoname, 1) == string.(s_autoname.u[])
@test_throws BoundsError ControlSystemsBase.output_names(s_autoname, 2)

s3 = NamedStateSpace(ssdata(G1)..., s1.x, s1.u, s1.y)
@test s3 == s1

@test s1[:y, :u] == s1
@test s1[[:y], [:u]] == s1

G3 = ControlSystemsBase.ssrand(1,2,3)
s3 = named_ss(G3, x = [:x1, :x2, :x3], u = [:u1, :u2], y=[:y])

@test s3[:y, :u1] == named_ss(G3[1,1], x = [:x1, :x2, :x3], u = [:u1], y=[:y])
@test s3[:y, :u2] == named_ss(G3[1,2], x = [:x1, :x2, :x3], u = [:u2], y=[:y])
@test s3[1, 2] == named_ss(G3[1,2], x = [:x1, :x2, :x3], u = [:u2], y=[:y]) # indexing with integer preserves names

@test s3[:y, [:u1, :u2]] == s3

G3 = partition(ssrand(2,3,4), 1, 1)
s3 = named_ss(G3, x = :x, u = :u, y=:y, z=:z, w=:w)
@test s3[:z, :w].sys == G3[1,1]
@test length(s3.u) == 3
@test length(s3.y) == 2

# Test prefix matching
G4 = ControlSystemsBase.ssrand(1,2,3)
s4 = named_ss(G4, x = :x, u = [:u1, :u2], y=[:y])
@test RobustAndOptimalControl.names2indices(:u, s4.u) == 1:2

G4 = ControlSystemsBase.ssrand(1,3,3)
s4 = named_ss(G4, x = [:x1, :x2, :x3], u = [:u1, :u2, :not_u], y=[:y])
@test RobustAndOptimalControl.names2indices(:u, s4.u) == 1:2

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

    G3 = 2*s2
    @test ss(G3) == 2*ss(s2)
    @test G3.u == s2.u
    @test G3.y != s2.y

    G3 = s2*2
    @test ss(G3) == ss(s2)*2
    @test G3.u != s2.u
    @test G3.y == s2.y

    G3 = s2/2.0
    @test ss(G3) == ss(s2)/2.0

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

    G = measure(s1, :x)
    @test G.C == ones(1, 1)
    @test G.y == [:x]
end

G1 = ss(1,1,1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
s2 = named_ss(G2, x = :z, u = :u2, y=:y2)
@test_nowarn @check_all_unique s1 s2

s1e = ExtendedStateSpace(s1, y=s1.y, u=s1.u, z=[], w=[])
@test s1.sys == system_mapping(s1e)
@test ss(s1e) == ss(s1)

s1d = c2d(s1, 1.0)
@test s1d.sys == c2d(s1.sys, 1.0)
@test s1d.x == s1.x
@test s1d.u == s1.u
@test s1d.y == s1.y
@test_nowarn plot(step(s1d))

fb = feedback(s1, s2)
# @test fb.x == [:x, :z]
@test fb.u == s1.u
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


# P = ssrand(1,2,3)
# addP1 = named_ss(ss([1  1], P.timeevol), u=[:vf_a, :yL], y=:uP)
# addL1 = named_ss(ss([I(P.nx) -I(P.nx)], P.timeevol), u=[:xr^P.nx; :xh^P.nx], y=:x_diff)
# addP = sumblock(:(uP = vf_a + yL))
# addL = sumblock(:(x_diff = xr - xh); n=P.nx)
# @test addP1 == addP
# @test addL1 == addL



# P = ssrand(1,2,3; Ts=1)
# addP1 = named_ss(ss([1  1], P.timeevol), u=[:vf_a, :yL], y=:uP)
# addL1 = named_ss(ss([I(P.nx) -I(P.nx)], P.timeevol), u=[:xr^P.nx; :xh^P.nx], y=:x_diff)
# addP = sumblock(:(uP = vf_a + yL); P.Ts)
# addL = sumblock(:(x_diff = xr - xh); n=P.nx, P.Ts)
# @test addP1 == addP
# @test addL1 == addL


P = ssrand(1,2,3)
addP1 = named_ss(ss([1  1], P.timeevol), u=[:vf_a, :yL], y=:uP)
addL1 = named_ss(ss([I(P.nx) -I(P.nx)], P.timeevol), u=[:xr^P.nx; :xh^P.nx], y=:x_diff)
addP = sumblock("uP = vf_a + yL")
addL = sumblock("x_diff = xr - xh"; n=P.nx)
@test addP1 == addP
@test addL1 == addL



P = ssrand(1,2,3; Ts=1)
addP1 = named_ss(ss([1  1], P.timeevol), u=[:vf_a, :yL], y=:uP)
addL1 = named_ss(ss([I(P.nx) -I(P.nx)], P.timeevol), u=[:xr^P.nx; :xh^P.nx], y=:x_diff)
addP = sumblock("uP = vf_a + yL"; P.Ts)
addL = sumblock("x_diff = xr - xh"; n=P.nx, P.Ts)
@test addP1 == addP
@test addL1 == addL

s4 = NamedStateSpace(ssdata(addP1.sys)..., addP1.Ts, addP1.x, addP1.u, addP1.y)
@test s4 == addP1

# Without spaces
addP = sumblock("uP=vf_a+yL"; P.Ts)
addL = sumblock("x_diff=xr-xh"; n=P.nx, P.Ts)
@test addP1 == addP
@test addL1 == addL

# More spaces
addP = sumblock("uP  =  vf_a  +  yL"; P.Ts)
addL = sumblock("x_diff  =  xr  -  xh"; n=P.nx, P.Ts)
@test addP1 == addP
@test addL1 == addL


# many
addmany = sumblock("uP  =  vf_a  +  yL + v"; P.Ts)
@test addmany.y == [:uP]
@test addmany.u == [:vf_a, :yL, :v]
@test addmany.D == [1 1 1]


# many
addmany = sumblock("uP  =  vf_a  +  yL - v"; P.Ts, n=2)
@test addmany.y == :uP^2
@test addmany.u == [:vf_a^2; :yL^2; :v^2]
@test addmany.D == [I(2) I(2) -I(2)]

##

G1 = ss(1,[1 1],1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x = [:x], u = :u, y=:y1) # test providing only symbol
s2 = named_ss(G2, x = [:z], u = [:u], y=[:y2])

connections = [:y1 => :u, :y2 => :u1]
w1 = [:u2]
G = connect([s1,s2], connections; w1)
@test G.u == w1
@test G.y == :y^2
@test G.sys == ss(ones(2,2), [1,0], I(2), 0)

sp = splitter(:u2, 2)
@test all(sp.D .== [1, 1])
@test sp.u == [:u2]
@test sp.y == :u2^2


@testset "complicated_feedback" begin
    @info "Testing complicated_feedback"
    include("../examples/complicated_feedback.jl")
end
