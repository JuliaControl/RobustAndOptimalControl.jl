module RobustAndOptimalControl

using LinearAlgebra, Statistics
using RecipesBase
using ControlSystems
import ControlSystems: ss, ssdata, ninputs, noutputs, nstates, isdiscrete, iscontinuous, to_matrix, timeevol, _string_mat_with_headers, PartionedStateSpace, common_timeevol

using ComponentArrays

using MonteCarloMeasurements, Optim, UnPack
import Distributions: Uniform

export ExtendedStateSpace
include("ExtendedStateSpace.jl")

export hinfsynthesize, hinfassumptions, hinfpartition, hinfsignals, bilinearc2d, bilineard2c, fudge_inv
include("hinfinity_design.jl")

include("plotting.jl")

export frequency_weighted_reduction, controller_reduction
include("reduction.jl")

export h2synthesize
include("h2_design.jl")

export LQGProblem, sensitivity, input_sensitivity, output_sensitivity, comp_sensitivity, input_comp_sensitivity, output_comp_sensitivity, controller, ff_controller, extended_controller, closedloop, static_gain_compensation, G_PS, G_CS, gangoffour, loopgain, stabilityrobustness, returndifference
include("lqg.jl")

export named_ss
include("named_systems.jl")

export find_lft, δ
include("find_lft.jl")


end
