module RobustAndOptimalControl

using ControlSystemsBase: issiso
using LinearAlgebra, Statistics
using RecipesBase
using UUIDs # To check if Plots is loaded in gofplot
using ControlSystemsBase
import ControlSystemsBase: ss, ssdata, ninputs, noutputs, nstates, isdiscrete, iscontinuous, to_matrix, timeevol, _string_mat_with_headers, common_timeevol, gangoffour, sensitivity,
input_sensitivity,
output_sensitivity,
comp_sensitivity,
input_comp_sensitivity,
output_comp_sensitivity,
G_PS,
G_CS

using ComponentArrays

using MonteCarloMeasurements, Optim, UnPack
import Distributions: Uniform
# import IntervalArithmetic
# import IntervalArithmetic: Interval
using GenericSchur
# using GenericLinearAlgebra # Causes issues with hessenberg due to pirate methods, the user must load this lib manually

import MatrixPencils, MatrixEquations, DescriptorSystems
using DescriptorSystems: dss

using ChainRulesCore

export show_construction, vec2sys
include("utils.jl")

export dss, hinfnorm2, linfnorm2, h2norm, hankelnorm, nugap, νgap, baltrunc2, stab_unstab, baltrunc_unstab, baltrunc_coprime
include("descriptor.jl")

export ExtendedStateSpace, system_mapping, performance_mapping, noise_mapping, ssdata_e, partition, ss
include("ExtendedStateSpace.jl")

export modal_form, schur_form, hess_form
include("canonical.jl")

export δ, δr, δc, δss, nominal, UncertainSS, uss, blocksort, sys_from_particles, ss2particles
include("uncertainty_interface.jl")

export makeweight, neglected_delay, gain_and_delay_uncertainty, neglected_lag, fit_complex_perturbations
include("weights.jl")

export hinfsynthesize, hinfassumptions, hinfpartition, hinfsignals, bilinearc2d, bilineard2c, fudge_inv, hinfgrad
include("hinfinity_design.jl")

export muplot, mvnyquistplot, specificationplot
include("plotting.jl")

export frequency_weighted_reduction, hsvd, controller_reduction, controller_reduction_weight, controller_reduction_plot
include("reduction.jl")

export h2synthesize
include("h2_design.jl")

export LQGProblem, sensitivity, input_sensitivity, output_sensitivity, comp_sensitivity, input_comp_sensitivity, output_comp_sensitivity, feedback_control, ff_controller, extended_controller, closedloop, static_gain_compensation, G_PS, G_CS, gangoffour
export lqr3, dare3
include("lqg.jl")

export NamedStateSpace, named_ss, expand_symbol, measure, connect, sumblock, splitter
include("named_systems2.jl")

export find_lft, δ
include("find_lft.jl")

export add_disturbance, add_low_frequency_disturbance, add_measurement_disturbance, add_resonant_disturbance, add_output_differentiator, add_output_integrator, add_input_integrator, add_input_differentiator
include("model_augmentation.jl")

export glover_mcfarlane, glover_mcfarlane_2dof, hanus, extended_gangoffour, ncfmargin
include("glover_mcfarlane.jl")


export diskmargin, Diskmargin, Disk, sim_diskmargin, loop_diskmargin, structured_singular_value, broken_feedback, robstab, loop_scaling, loop_scale, passivity_index, ispassive, passivityplot
include("diskmargin.jl")
include("mimo_diskmargin.jl")


export nu_reduction, nu_reduction_recursive
include("mcm_nugap.jl")

end
