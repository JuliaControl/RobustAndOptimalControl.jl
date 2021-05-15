using RobustAndOptimalControl
using Test

@testset "RobustAndOptimalControl.jl" begin
    # Write your tests here.
    @testset "H∞ design" begin
        @info "Testing hinf_design"
        include("test_hinf_design.jl")
    end

    @testset "H2 design" begin
        @info "Testing H2 design"
        include("test_h2_design.jl")
    end

    @testset "Named systems" begin
        @info "Testing Named systems"
        include("test_named_systems.jl")
    end
end
