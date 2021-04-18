module RobustAndOptimalControl

using LinearAlgebra, Statistics
using RecipesBase
using ControlSystems
import ControlSystems: ss, ssdata, ninputs, noutputs, nstates, isdiscrete, iscontinuous

export hinfsynthesize, hinfassumptions, hinfpartition, hinfsignals, bilinearc2d, bilineard2c

include("ExtendedStateSpace.jl")
include("hinfinity_design.jl")
include("plotting.jl")


end
