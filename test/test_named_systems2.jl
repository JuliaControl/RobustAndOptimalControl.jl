using RobustAndOptimalControl, ControlSystemsBase, LinearAlgebra, Test, Plots
using RobustAndOptimalControl: check_unique, check_all_unique

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

@test_nowarn RobustAndOptimalControl.show_construction(s_autoname)

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

        @test op(I(1), s1).sys == op(I(1), G1)
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

    s11 = s1*s1
    @test allunique(s11.x)

    @test ControlSystemsBase.numeric_type(big(1.0)*s1) <: BigFloat
    @test ControlSystemsBase.numeric_type(s1*big(1.0)) <: BigFloat
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

    @test add_output(s1, [0.2]) isa NamedStateSpace
    @test add_output(s1, [0.2], y=[:hej]).y[2] === :hej

end

G1 = ss(1,1,1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
s2 = named_ss(G2, x = :z, u = :u2, y=:y2)
@test_nowarn check_all_unique(s1, s2)

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
@test_throws ArgumentError check_all_unique(s1, s2)

# Same x names, test automatic generation of new names
s1 = named_ss(ssrand(1,1,2), x = :x, u = :u1, y=:y1)
s2 = named_ss(ssrand(1,1,2), x = :x, u = :u2, y=:y2)
s12 = [s1; s2]
@test occursin("x1", string(s12.x[1]))
@test occursin("x2", string(s12.x[2]))
@test occursin("x1", string(s12.x[3]))
@test occursin("x2", string(s12.x[4]))

# When systems have names, use these to create the new x names
s1 = named_ss(ssrand(1,1,2), "P", x = :x, u = :u1, y=:y1)
s2 = named_ss(ssrand(1,1,2), "C", x = :x, u = :u2, y=:y2)
s12 = [s1; s2]
@test s12.x[1] == :Px1
@test s12.x[2] == :Px2
@test s12.x[3] == :Cx1
@test s12.x[4] == :Cx2

# When one system is missing a name, we use the existing name only
s1 = named_ss(ssrand(1,1,2), "P", x = :x, u = :u1, y=:y1)
s2 = named_ss(ssrand(1,1,2), x = :x, u = :u2, y=:y2)
s12 = [s1; s2]
@test s12.x[1] == :Px1
@test s12.x[2] == :Px2
@test s12.x[3] == :x1
@test s12.x[4] == :x2

# If more than one system is missing a name, we do not use the system names
s1 = named_ss(ssrand(1,1,2), "P", x = :x, u = :u1, y=:y1)
s2 = named_ss(ssrand(1,1,2), x = :x, u = :u2, y=:y2)
s3 = named_ss(ssrand(1,1,2), x = :x, u = :u3, y=:y3)
s12 = [s1; s2; s3]
@test occursin("x1", string(s12.x[1]))
@test occursin("x2", string(s12.x[2]))
@test occursin("x1", string(s12.x[3]))
@test occursin("x2", string(s12.x[4]))
@test occursin("x1", string(s12.x[5]))
@test occursin("x2", string(s12.x[6]))


op = RobustAndOptimalControl.operating_point(s1)
@test op.x == zeros(s1.nx)
@test op.u == zeros(s1.nu)

opx = randn(s1.nx)
opu = randn(s1.nu)
RobustAndOptimalControl.set_operating_point!(s1, (x = opx, u = opu))
@test RobustAndOptimalControl.operating_point(s1) == (x = opx, u = opu)

s1b, T = balance_statespace(s1)
opb = RobustAndOptimalControl.infer_operating_point(s1b, s1)

@test RobustAndOptimalControl.operating_point(s1b).x ≈ opb.x

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

    G1 = named_ss(ssrand(1,1,1, Ts=1), "G1")
    G2 = named_ss(ssrand(1,1,1, Ts=1), "G2")
    gangoffourplot(G1, G2) # tests some convert methods for I to discrete


    G1 = named_ss(ssrand(1,1,1), "G1")
    # Scalars
    @test_nowarn G1*1
    @test_nowarn 1*G1

    # Transfer function
    @test (G1*tf(1, [1,1])).sys == (G1*ss(tf(1, [1,1]))).sys
    @test (tf(1, [1,1])*G1).sys == (ss(tf(1, [1,1]))*G1).sys

    # Matrix
    @test (G1*ones(1,1)).sys == (G1*ss(ones(1,1))).sys
    @test (ones(1,1)*G1).sys == (ss(ones(1,1))*G1).sys

    # if the matrix is diagonal, the names are `u_scaled`
    @test endswith(string((G1*ones(1,1)).u[]), "_scaled")

    # If the matrix is not diagonal, the names are generic
    G1 = named_ss(ssrand(1,2,1), "G1")
    @test !endswith(string((G1*ones(2,2)).u[1]), "_scaled")
end


## Test signal names after feedback etc.

P = named_ss(ssrand(1,1,1), "P")
C = named_ss(ssrand(1,1,1), "C")
S, PS, CS, T = gangoffour(P, C)

@test PS.u[] == :Pu
@test PS.y[] == :Py
@test CS.u[] == :Cu
@test CS.y[] == :Cy
@test T.u[] == :Cu
@test T.y[] == :Py
@test S.u[] == :Py
@test S.y[] == :Cu

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
addP1 = named_ss(ss([1  1], P.timeevol), u=[:vf_a, :yL], y=:uP, name="sumblock")
addL1 = named_ss(ss([I(P.nx) -I(P.nx)], P.timeevol), u=[:xr^P.nx; :xh^P.nx], y=:x_diff, name="sumblock")
addP = sumblock("uP = vf_a + yL")
addL = sumblock("x_diff = xr - xh"; n=P.nx)
@test addP1 == addP
@test addL1 == addL



P = ssrand(1,2,3; Ts=1)
addP1 = named_ss(ss([1  1], P.timeevol), u=[:vf_a, :yL], y=:uP, name="sumblock")
addL1 = named_ss(ss([I(P.nx) -I(P.nx)], P.timeevol), u=[:xr^P.nx; :xh^P.nx], y=:x_diff, name="sumblock")
addP = sumblock("uP = vf_a + yL"; P.Ts)
addL = sumblock("x_diff = xr - xh"; n=P.nx, P.Ts)
@test addP1 == addP
@test addL1 == addL

s4 = NamedStateSpace(ssdata(addP1.sys)..., addP1.Ts, addP1.x, addP1.u, addP1.y, addP1.name)
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

## Test bug with prefix matching, where one of the provided indices matches more than one variable (:y), but there are still remaining indices to find (:za)
G1 = ssrand(4,2,2)
s1 = named_ss(G1, u = :u, y=[:y1, :y2, :w, :za]) 
@test s1[[:y, :z], :] == s1[[1,2,4], :]

## Test conversion of tf
nss = named_ss(tf(1, [1,1]))
@test nss isa NamedStateSpace



## Test where one external input goes to several inputs with the same name, with and without using a splitter
s1 = ssrand(1,2,2)
s2 = ssrand(1,2,2)
sys_1=named_ss(s1, u=[:in_x, :in_y], y=:out_z, x=[:x1, :x2])
sys_2=named_ss(s2, u=[:in_y, :in_z], y=:out, x = [:x3, :x4])

@test_throws "u names not unique. Repeated names: [:in_y] To allow connecting a single input signal to several inputs with the same name, pass `unique = false`." connect([sys_1, sys_2], [:out_z => :in_z]; w1 = [:in_x, :in_y], z1 = :out)
sys_connect = connect([sys_1, sys_2], [:out_z => :in_z]; w1 = [:in_x, :in_y], z1 = :out, unique=false)

sys_1=named_ss(s1, u=[:in_x, :in_y1], y=:out_z, x=[:x1, :x2])
sys_2=named_ss(s2, u=[:in_y2, :in_z], y=:out, x = [:x3, :x4])
split = splitter(:in_y, 2)
sys_connect2 = connect([sys_1, sys_2, split], [:out_z => :in_z, :in_y1=>:in_y1, :in_y2=>:in_y2]; w1 = [:in_x, :in_y], z1 = :out)

# Make gensym state names equal to assist equality testing below
sys_connect2.x .= sys_connect.x
@test sys_connect ≈ sys_connect2

# Discrete
s1 = ssrand(1,2,2, Ts=1)
s2 = ssrand(1,2,2, Ts=1)
sys_1=named_ss(s1, u=[:in_x, :in_y], y=:out_z, x=[:x1, :x2])
sys_2=named_ss(s2, u=[:in_y, :in_z], y=:out, x = [:x3, :x4])

@test_throws "u names not unique. Repeated names: [:in_y] To allow connecting a single input signal to several inputs with the same name, pass `unique = false`." connect([sys_1, sys_2], [:out_z => :in_z]; w1 = [:in_x, :in_y], z1 = :out)
sys_connect = connect(
    [sys_1, sys_2],
    [:out_z => :in_z];
    w1 = [:in_x, :in_y],
    z1 = :out,
    unique=false
)

sys_1=named_ss(s1, u=[:in_x, :in_y1], y=:out_z, x=[:x1, :x2])
sys_2=named_ss(s2, u=[:in_y2, :in_z], y=:out, x = [:x3, :x4])
split = splitter(:in_y, 2, sys_1.timeevol)
sys_connect2 = connect([sys_1, sys_2, split], [:out_z => :in_z, :in_y1=>:in_y1, :in_y2=>:in_y2]; w1 = [:in_x, :in_y], z1 = :out)

# Make gensym state names equal to assist equality testing below
sys_connect2.x .= sys_connect.x
@test sys_connect ≈ sys_connect2


## Multiple identical names
s1 = ssrand(1,3,2)
s2 = ssrand(1,3,2)
s3 = ssrand(1,1,1)

sys_1 = named_ss(s1, u=[:in_a, :in_b, :in_c], y=[:ui_ref], x=[:x1, :x2])
sys_2 = named_ss(s2, u=[:in_e, :in_c, :in_d], y=:out_a, x=[:x3, :x4]);
sys_3 = named_ss(s3, u=:in_d, y=:out_e, x=:x5);
sys_connect = connect([sys_1, sys_2, sys_3],
    [:out_a => :in_a];
    w1 = [:in_b, :in_c, :in_d],
    z1 = :ui_ref,
    unique = false
)


split1 = splitter(:in_c, 2, sys_1.timeevol)
split2 = splitter(:in_d, 2, sys_1.timeevol)

sys_1 = named_ss(s1, u=[:in_a, :in_b, :in_c1], y=[:ui_ref], x=[:x1, :x2])
sys_2 = named_ss(s2, u=[:in_e, :in_c2, :in_d1], y=:out_a, x=[:x3, :x4]);
sys_3 = named_ss(s3, u=:in_d2, y=:out_e, x=:x5);
sys_connect2 = connect([sys_1, sys_2, sys_3, split1, split2],
    [:out_a => :in_a, :in_c1=>:in_c1, :in_c2=>:in_c2, :in_d1=>:in_d1, :in_d2=>:in_d2];
    w1 = [:in_b, :in_c, :in_d],
    z1 = :ui_ref,
    unique = false
)
sys_connect2.x .= sys_connect.x
@test sys_connect ≈ sys_connect2


##
s1 = ssrand(1,2,2)
s2 = ssrand(1,2,2)
sys1 = named_ss(s1, "A", u=[:u1, :u2])
sys2 = named_ss(s2, "B", u=[:u2, :u3])

sys12 = connect(
    [sys1, sys2],
    Pair{Symbol, Symbol}[];
    w1 = [:u1, :u2, :u3],
    z1 = [sys1.y; sys2.y],
    unique=false
)

@test sys12.B[1:2,1] == sys1.B[:, 1]
@test sys12.B[:,2] == [sys1.B[:, 2]; sys2.B[:, 1]]
@test sys12.B[3:4,3] == sys2.B[:, 2]

## Inv
s1 = named_ss(ssrand(2,2,2))
isys = inv(s1)
@test isys.sys == inv(s1.sys)
@test isys.x == s1.x
@test isys.u == s1.y # Names are reversed
@test isys.y == s1.u

isys = 2/s1
@test isys.sys == 2/s1.sys
@test isys.x != s1.x