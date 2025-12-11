if haskey(ENV, "CI")
    ENV["PLOTS_TEST"] = "true"
    ENV["GKSwstype"] = "100" # gr segfault workaround
end
using Plots
using RobustAndOptimalControl
using LinearAlgebra
using Test


@testset "RobustAndOptimalControl.jl" begin

    @testset "extendedstatespace" begin
        @info "Testing extendedstatespace"
        include("test_extendedstatespace.jl")
    end

    @testset "utils" begin
        @info "Testing utils"
        include("test_utils.jl")
    end

    @testset "descriptor" begin
        @info "Testing descriptor"
        include("test_descriptor.jl")
    end

    @testset "uncertainty" begin
        @info "Testing uncertainty"
        include("test_uncertainty.jl")
    end

    @testset "diskmargin" begin
        @info "Testing diskmargin"
        include("test_diskmargin.jl")
    end

    @testset "weights" begin
        @info "Testing weights"
        include("test_weights.jl")
    end

    @testset "hinfpartition" begin
        @info "Testing hinfpartition"
        include("test_hinfpartition.jl")
    end

    @testset "H∞ design" begin
        @info "Testing hinf_design"
        include("test_hinf_design.jl")
    end

    @testset "manual_hinf" begin
        @info "Testing manual_hinf"
        include("test_manual_hinf.jl")
    end

    @testset "H2 design" begin
        @info "Testing H2 design"
        include("test_h2_design.jl")
    end

    @testset "LQG" begin
        @info "Testing LQG"
        include("test_lqg.jl")
    end

    @testset "LQI" begin
        @info "Testing LQI"
        include("test_lqi.jl")
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

    @testset "glover_mcfarlane" begin
        @info "Testing glover_mcfarlane"
        include("test_glover_mcfarlane.jl")
    end

    @testset "hinfgrad" begin
        @info "Testing hinfgrad"
        include("test_hinfgrad.jl")
    end

    # The following test is run last. It loads GenericLinearAlgebra which has pirate linear algebra methods. Those methods screw up other tests.
    @testset "high_precision_hinf" begin
        @info "Testing high_precision_hinf"
        include("test_high_precision_hinf.jl")
    end

    @testset "canonical_forms" begin
        @info "Testing canonical_forms"
        include("test_canonical_forms.jl")
    end

    @testset "complicated feedback exmaple" begin
        @info "Testing complicated feedback exmaple"
        include("../examples/complicated_feedback.jl")
    end

    @testset "Flexible servo" begin
        @info "Testing Flexible servo"
        include("../examples/flexible_servo.jl")
    end

    @testset "Uncertain" begin
        @info "Testing Flexible servo"
        include("../examples/uncertain.jl")
    end

    @testset "lqg_mpc_disturbance" begin
        @info "Testing LQG MPC example"
        include("../examples/lqg_mpc_disturbance.jl")
    end

    @testset "mcm_nugap" begin
        @info "Testing mcm_nugap"
        include("test_mcm_nugap.jl")
    end

    @testset "LinearMPC extension" begin
        @info "Testing LinearMPC extension"
        include("test_linearmpc_ext.jl")
    end

end
