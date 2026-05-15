using Test
using ControlSystemsBase
using RobustAndOptimalControl
using LinearAlgebra
using Plots


G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
Q = diagm([0,5])
R = [1.0;;]
K = kalman(G,Q,R)
obs = observer_predictor(G,K; output_state=true)

Q1 = diagm([0.488,0,100])
Q2 = [1/100;;]
L = lqi(G,Q1,Q2)

C0 = RobustAndOptimalControl.lqi_controller(G, obs, Q1, Q2)

@test C0.nu == 2
@test 0 ∈ poles(C0)

Gn = named_ss(G)
H = feedback(C0, Gn, w1 = :y_plant_r, z2=Gn.y, u1=:y_plant, pos_feedback=true)

@test dcgain(H)[2] ≈ 1

res = step(H, 50)
@test res.y[:, end] ≈ [dcgain(feedback(-C0[:, 2], G))[]; 1.0]
# plot(res)


# MIMO continuous, all outputs integrated
G2 = ss([-1.0 0.2; 0.1 -2.0], [1.0 0; 0 1.0], [1.0 0; 0 1.0], 0)
K2 = kalman(G2, I(G2.nx), I(G2.ny))
obs2 = observer_predictor(G2, K2; output_state=true)
Q1_mimo = diagm([1.0, 1.0, 10.0, 10.0])
Q2_mimo = Matrix{Float64}(I(G2.nu))
C_mimo = RobustAndOptimalControl.lqi_controller(G2, obs2, Q1_mimo, Q2_mimo)
@test C_mimo.nu == 2*G2.ny       # [r; y]
@test C_mimo.ny == G2.nu
# Each output should have an integral mode in the controller
@test count(p -> abs(p) < 1e-6, poles(C_mimo)) == G2.ny

ref_syms_mimo = Symbol.("y_plant" .* string.(1:G2.ny) .* "_r")
y_plant_syms = Symbol.("y_plant" .* string.(1:G2.ny))
G2n = named_ss(G2, name="plant", x=:x_plant, y=:y_plant, u=:u_plant)
H2 = feedback(C_mimo, G2n, w1 = ref_syms_mimo,
              z2 = G2n.y, u1 = y_plant_syms, pos_feedback = true)
@test isstable(minreal(H2))
# Reference-to-output DC gain should be identity on the plant outputs.
# H2 outputs the plant outputs (z2 = G2n.y) and any external outputs of C_mimo (its u).
H2_dc = dcgain(H2)
# Extract the y-subset of outputs by name
y_out_inds = [findfirst(==(s), H2.y) for s in G2n.y]
@test H2_dc[y_out_inds, :] ≈ I(G2.ny) atol=1e-8


# Partial integrator_outputs: integrate only output 1 on a 2-output plant
Q1_partial = diagm([1.0, 1.0, 10.0])  # nx + 1 integrator
C_partial = RobustAndOptimalControl.lqi_controller(G2, obs2, Q1_partial, Q2_mimo; integrator_outputs=[1])
# Controller takes [r_for_integrated_output; all y] = [1 ref; ny measurements]
@test C_partial.nu == 1 + G2.ny
@test count(p -> abs(p) < 1e-6, poles(C_partial)) == 1
H_partial = feedback(C_partial, G2n, w1 = [ref_syms_mimo[1]],
                     z2 = G2n.y, u1 = y_plant_syms, pos_feedback = true)
@test isstable(minreal(H_partial))
H_partial_dc = dcgain(H_partial)
y_out_inds_p = [findfirst(==(s), H_partial.y) for s in G2n.y]
# Only the integrated output should track exactly
@test H_partial_dc[y_out_inds_p[1], 1] ≈ 1 atol=1e-8


# Discrete SISO with explicit ϵ
Ts = 0.1
Gd = c2d(G, Ts)
Kd = kalman(Gd, Q, R)
obsd = observer_predictor(Gd, Kd; output_state=true)
Q1d = diagm([0.488, 0, 100.0])
Q2d = [1/100;;]
ϵd = 1e-4
Cd = RobustAndOptimalControl.lqi_controller(Gd, obsd, Q1d, Q2d; ϵ=ϵd)
@test Cd.nu == 2
# Discrete integrator pole near 1 (offset by ϵ)
@test any(p -> abs(p - (1 - ϵd)) < 1e-6, poles(Cd))
Gdn = named_ss(Gd)
Hd = feedback(Cd, Gdn, w1 = :y_plant_r, z2=Gdn.y, u1=:y_plant, pos_feedback=true)
@test isstable(minreal(Hd))
@test dcgain(Hd)[2] ≈ 1 atol=1e-2


# LQGProblem method
prob = LQGProblem(G, diagm([0.488, 0]), [1/100;;], Matrix{Float64}(Q), Matrix{Float64}(R))
Qi = [100.0;;]
C_prob = RobustAndOptimalControl.lqi_controller(prob, Qi)
@test C_prob.nu == 2
@test 0 ∈ poles(C_prob)
# Dimensions should match the explicit-observer form
@test C_prob.nx == C0.nx
@test C_prob.ny == C0.ny
