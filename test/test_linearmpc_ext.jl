using Test
using LinearAlgebra
using LinearMPC
using RobustAndOptimalControl

# Test basic LQGProblem conversion (from docstring example)
Ts = 0.1
sys = ss([1 Ts; 0 1], [0; Ts], [1 0], 0, Ts)

# Create LQGProblem: Q1=state cost, Q2=input cost, R1=process noise, R2=measurement noise
lqg = LQGProblem(sys, I(2), I(1), I(2), 0.01*I(1))

# Convert to LinearMPC with constraints
mpc = LinearMPC.MPC(lqg; N=20, umin=[-0.3], umax=[0.3])

# Check basic properties
@test mpc.model.Ts == Ts
@test mpc.model.F == sys.A
@test mpc.model.G == sys.B

# Simulate and verify it works
sim = LinearMPC.Simulation(mpc; N=100, r=[1.0, 0])


# Check constraint satisfaction
@test all(sim.us .>= -0.3 - 1e-6)
@test all(sim.us .<= 0.3 + 1e-6)

# Check that output converges toward reference
final_output = (sys.C * sim.xs[:, end])[1]
@test abs(final_output - 1.0) < 0.1  # Should be close to reference

@testset "LQGProblem with Q3 (input rate penalty)" begin
    # Test with input rate cost
    Q3 = 0.5 * I(1)
    mpc_q3 = LinearMPC.MPC(lqg; N=15, Q3, umin=[-0.5], umax=[0.5])

    sim_q3 = LinearMPC.Simulation(mpc_q3; N=50, r=[0.5, 0])

    # With Q3, control should be smoother (smaller rate of change)
    du = diff(sim_q3.us, dims=2)
    max_rate = maximum(abs.(du))
    @test max_rate < 0.1  # Rate should be limited due to Q3 penalty
end

@testset "LQGProblem with terminal cost Qf" begin
    # Test with terminal cost
    Qf = 10.0 * I(2)  # Higher terminal cost
    mpc_qf = LinearMPC.MPC(lqg; N=10, Qf, umin=[-1.0], umax=[1.0])

    sim_qf = LinearMPC.Simulation(mpc_qf; N=30, r=[0.3, 0])

    # Should still work and converge
    final_output = (sys.C * sim_qf.xs[:, end])[1]
    @test abs(final_output - 0.3) < 0.1
end

@testset "LQGProblem MIMO system" begin
    # Test with 2-input 2-output system
    A = [0.9 0.1; 0.05 0.95]
    B = [1.0 0.0; 0.0 1.0]
    C = [1.0 0.0; 0.0 1.0]
    D = zeros(2, 2)
    sys_mimo = ss(A, B, C, D, Ts)

    lqg_mimo = LQGProblem(sys_mimo, I(2), I(2), I(2), 0.01*I(2))
    mpc_mimo = LinearMPC.MPC(lqg_mimo; N=15,
                                umin=[-0.5, -0.5], umax=[0.5, 0.5])

    sim_mimo = LinearMPC.Simulation(mpc_mimo; N=50, r=[0.3, -0.2])

    # Check dimensions
    @test size(sim_mimo.us, 1) == 2
    @test size(sim_mimo.xs, 1) == 2

    # Check constraints on both inputs
    @test all(sim_mimo.us[1, :] .>= -0.5 - 1e-6)
    @test all(sim_mimo.us[1, :] .<= 0.5 + 1e-6)
    @test all(sim_mimo.us[2, :] .>= -0.5 - 1e-6)
    @test all(sim_mimo.us[2, :] .<= 0.5 + 1e-6)
end

@testset "LQGProblem continuous-time error" begin
    # Test that continuous-time systems throw an error
    sys_cont = ss([0 1; -1 -1], [0; 1], [1 0], 0)  # Continuous-time
    lqg_cont = LQGProblem(sys_cont, I(2), I(1), I(2), 0.01*I(1))

    @test_throws ErrorException LinearMPC.MPC(lqg_cont; N=10)
end
