using RobustAndOptimalControl
using Test

@testset "RobustAndOptimalControl.jl" begin
    @testset "H∞ design" begin
        @info "Testing hinf_design"
        include("test_hinf_design.jl")
    end

    @testset "H2 design" begin
        @info "Testing H2 design"
        include("test_h2_design.jl")
    end

    @testset "LQG" begin
        @info "Testing LQG"
        include("test_lqg.jl")
    end

    @testset "Named systems" begin
        @info "Testing Named systems"
        include("test_named_systems.jl")
    end

    @testset "find_lft" begin
        @info "Testing find_lft"
        include("test_find_lft.jl")
    end
end
