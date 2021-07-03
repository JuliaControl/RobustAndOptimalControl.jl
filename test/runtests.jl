using RobustAndOptimalControl
using Test

@testset "RobustAndOptimalControl.jl" begin

    @testset "hinfpartition" begin
        @info "Testing hinfpartition"
        include("test_hinfpartition.jl")
    end

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
        include("test_named_systems2.jl")
    end

    @testset "find_lft" begin
        @info "Testing find_lft"
        include("test_find_lft.jl")
    end

    @testset "reduction" begin
        @info "Testing test_reduction"
        include("test_reduction.jl")
    end

    @testset "augmentation" begin
        @info "Testing augmentation"
        include("test_augmentation.jl")
    end
end
