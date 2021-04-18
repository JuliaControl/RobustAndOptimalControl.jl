using RobustAndOptimalControl
using Test

@testset "RobustAndOptimalControl.jl" begin
    # Write your tests here.
    @testset "H∞ design" begin
        @info "Testing hinf_design"
        include("test_hinf_design.jl")
    end
end
