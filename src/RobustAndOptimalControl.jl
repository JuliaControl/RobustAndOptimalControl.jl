module RobustAndOptimalControl

using LinearAlgebra, Statistics
using RecipesBase
using ControlSystems
import ControlSystems: ss, ssdata, ninputs, noutputs, nstates, isdiscrete, iscontinuous, to_matrix, timeevol, _string_mat_with_headers, PartionedStateSpace, common_timeevol

export ExtendedStateSpace

export hinfsynthesize, hinfassumptions, hinfpartition, hinfsignals, bilinearc2d, bilineard2c, fudge_inv

export h2synthesize

export frequency_weighted_reduction, controller_reduction

export named_ss

include("ExtendedStateSpace.jl")
include("hinfinity_design.jl")
include("plotting.jl")
include("reduction.jl")
include("h2_design.jl")
include("named_systems.jl")


end
