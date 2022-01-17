module RobustAndOptimalControl

using ControlSystems: issiso
using LinearAlgebra, Statistics
using RecipesBase
using ControlSystems
import ControlSystems: ss, ssdata, ninputs, noutputs, nstates, isdiscrete, iscontinuous, to_matrix, timeevol, _string_mat_with_headers, PartionedStateSpace, common_timeevol

using ComponentArrays

using MonteCarloMeasurements, Optim, UnPack
import Distributions: Uniform
import IntervalArithmetic
import IntervalArithmetic: Interval


import MatrixPencils, MatrixEquations, DescriptorSystems
using DescriptorSystems: dss

using ChainRulesCore

export show_construction, vec2sys
include("utils.jl")

export dss, hinfnorm2, h2norm, hankelnorm, nugap, νgap, baltrunc2
include("descriptor.jl")

export ExtendedStateSpace, system_mapping, performance_mapping, noise_mapping, ssdata_e, partition, ss
include("ExtendedStateSpace.jl")

export δ, δr, δc, δss, nominal, UncertainSS, uss, blocksort
include("uncertainty_interface.jl")

export makeweight, neglected_delay, gain_and_delay_uncertainty, neglected_lag, fit_complex_perturbations
include("weights.jl")

export hinfsynthesize, hinfassumptions, hinfpartition, hinfsignals, bilinearc2d, bilineard2c, fudge_inv, hinfgrad
include("hinfinity_design.jl")

export muplot, mvnyquistplot, specificationplot
include("plotting.jl")

export frequency_weighted_reduction, hsvd
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

export glover_mcfarlane, glover_mcfarlane_2dof, hanus, extended_gangoffour, ncfmargin
include("glover_mcfarlane.jl")


export diskmargin, Diskmargin, Disk, sim_diskmargin, loop_diskmargin, structured_singular_value, broken_feedback, robstab, loop_scaling, loop_scale
include("diskmargin.jl")
include("mimo_diskmargin.jl")

end
