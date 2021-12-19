module RobustAndOptimalControl

using ControlSystems: issiso
using LinearAlgebra, Statistics
using RecipesBase
using ControlSystems
import ControlSystems: ss, ssdata, ninputs, noutputs, nstates, isdiscrete, iscontinuous, to_matrix, timeevol, _string_mat_with_headers, PartionedStateSpace, common_timeevol

using ComponentArrays

using MonteCarloMeasurements, Optim, UnPack
import Distributions: Uniform

import MatrixPencils, MatrixEquations

export ExtendedStateSpace, system_mapping, performance_mapping, ssdata_e
include("ExtendedStateSpace.jl")

export Weights, makeweight
include("weights.jl")

export hinfsynthesize, hinfassumptions, hinfpartition, hinfsignals, bilinearc2d, bilineard2c, fudge_inv
include("hinfinity_design.jl")

include("plotting.jl")

export frequency_weighted_reduction, controller_reduction, hsvd
include("reduction.jl")

export h2synthesize
include("h2_design.jl")

export LQGProblem, sensitivity, input_sensitivity, output_sensitivity, comp_sensitivity, input_comp_sensitivity, output_comp_sensitivity, controller, ff_controller, extended_controller, closedloop, static_gain_compensation, G_PS, G_CS, gangoffour, loopgain, stabilityrobustness, returndifference
include("lqg.jl")

export NamedStateSpace, named_ss, expand_symbol, measure, connect, sumblock, splitter
include("named_systems2.jl")

export find_lft, δ
include("find_lft.jl")

export add_disturbance, add_low_frequency_disturbance, add_measurement_disturbance, add_resonant_disturbance, add_output_differentiator, add_output_integrator, add_input_integrator, add_input_differentiator
include("model_augmentation.jl")


end
