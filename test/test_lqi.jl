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
