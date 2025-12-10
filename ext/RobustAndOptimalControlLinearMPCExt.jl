"""
    RobustAndOptimalControlLinearMPCExt

This extension allows you to convert an `LQGProblem` from RobustAndOptimalControl.jl to a `LinearMPC.MPC` object from LinearMPC.jl.
"""
module RobustAndOptimalControlLinearMPCExt

using LinearMPC
using LinearAlgebra
using RobustAndOptimalControl: LQGProblem


"""
    LinearMPC.MPC(prob::LQGProblem; N, Nc=N, Q3=nothing, Qf=nothing,
                 umin=nothing, umax=nothing, ymin=nothing, ymax=nothing, kwargs...)

Convert an `LQGProblem` from RobustAndOptimalControl.jl to a `LinearMPC.MPC` object.

# Arguments
- `prob::LQGProblem`: The LQG problem to convert
- `N::Int`: Prediction horizon (required)
- `Nc::Int`: Control horizon (default: `N`)
- `Q3`: Input rate cost matrix (default: zeros)
- `Qf`: Terminal state cost matrix (default: none)
- `umin`, `umax`: Input bounds (default: none)
- `ymin`, `ymax`: Output bounds (default: none)
- `kwargs...`: Additional arguments passed to LinearMPC.MPC

# Notes
- Only discrete-time systems are supported
- The Kalman filter/observer from LQGProblem is not converted
- Uses C1 (performance output) as the controlled output matrix
- Cost matrices Q1, Q2 from LQGProblem map to Q, R in LinearMPC

# Example
```julia
using RobustAndOptimalControl, LinearMPC, LinearAlgebra

sys = ss([1 0.1; 0 1], [0; 0.1], [1 0], 0, 0.1)
lqg = LQGProblem(sys, I(2), I(1), I(2), 0.01*I(1))
mpc = LinearMPC.MPC(lqg; N=20, umin=[-0.3], umax=[0.3])

sim = LinearMPC.Simulation(mpc, N=100, r=[1.0, 0])

using Plots
plot(sim)
```
"""
function LinearMPC.MPC(prob::LQGProblem;
    N::Int,
    Nc::Int = N,
    Q3 = zeros(0,0),
    Qf = zeros(0,0),
    umin = zeros(0),
    umax = zeros(0),
    ymin = zeros(0),
    ymax = zeros(0),
    kwargs...
)
    # Validate discrete-time system
    if !ControlSystemsBase.isdiscrete(prob)
        error("Only discrete-time systems are supported. Got a continuous-time system.")
    end

    # Extract system matrices
    F = Matrix(prob.A)   # State transition matrix
    G = Matrix(prob.B2)  # Control input matrix
    C = Matrix(prob.C1)  # Performance output matrix

    # Get sampling time
    Ts = prob.Ts

    # Get dimensions
    nu = size(G, 2)
    ny = size(C, 1)

    # Create the LinearMPC.MPC object
    mpc = LinearMPC.MPC(F, G; Ts, C, Np=N, Nc, kwargs...)

    # Set objective
    Q = Matrix(prob.Q1)
    R = Matrix(prob.Q2)
    Rr = Matrix(Q3)

    LinearMPC.set_objective!(mpc; Q, R, Rr, Qf)
    LinearMPC.set_bounds!(mpc; umin, umax, ymin, ymax)

    # Set labels from system names
    # sys = prob.sys
    # set_labels!(mpc; x=Symbol.(state_names(sys)), u=Symbol.(input_names(sys)))

    return mpc
end

end # module