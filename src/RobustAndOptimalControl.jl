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

export named_ss
include("named_systems.jl")

export find_lft, δ
include("find_lft.jl")


end
