#=
This file demonstrates
- Adding integral action to an LQG controller by means of state augmentation
- Conversion of the LQG controller into an MPC contorller with added constraints
=#
using RobustAndOptimalControl, ControlSystemsBase, Plots, LinearAlgebra
Ts = 1 # Sample time
G = c2d(ss(tf(1, [10, 1])), Ts) # Process model

nx  = G.nx
nu  = G.nu
ny  = G.ny
x0  = zeros(G.nx) # Initial condition

Q1 = 100diagm(ones(G.nx)) # state cost matrix
Q2 = 0.01diagm(ones(nu))  # control cost matrix

R1 = 0.001I(nx) # State noise covariance
R2 = I(ny)      # measurement noise covariance
prob = LQGProblem(G, Q1, Q2, R1, R2)

disturbance = (x, t) -> t * Ts ≥ 10 # This is our load disturbance, a step at ``t = 10``
Gcl = G_PS(prob)  # This forms the transfer function from load disturbance to output
res = lsim(Gcl, disturbance, 100)
plot(res)

Gd1 = add_low_frequency_disturbance(G, ϵ = 1e-6, measurement=false)
Gd2 = add_low_frequency_disturbance(G, ϵ = 1e-6, measurement=true) # The ϵ moves the integrator pole slightly into the stable region
plots = map([Gd1, Gd2]) do Gd
    nx = Gd.nx

    C  = Gd.C
    Q1 = 100C'diagm(ones(G.nx)) * C # state cost matrix
    x0 = zeros(nx)

    R1 = diagm([0.001, 1])
    R2 = I(ny)
    prob = LQGProblem(Gd, Q1, Q2, R1, R2)
    Gcl  = [G_PS(prob); -comp_sensitivity(prob)] # -comp_sensitivity(prob) is the same as the transfer function from load disturbance to control signal
    res  = lsim(Gcl, disturbance, 100)
    f1 = plot(res, ylabel=["y" "u"]); ylims!((-0.05, 0.3), sp = 1)

    w = exp10.(LinRange(-3, log10(pi / Ts), 200))
    f2 = gangoffourplot(prob, w, lab = "", legend = :bottomright)

    R1 = diagm([0.001, 0.2]) # Reduce the noise on the integrator state from 1 to 0.2
    R2 = I(ny)
    prob = LQGProblem(Gd, Q1, Q2, R1, R2)

    Gcl = [G_PS(prob); -comp_sensitivity(prob)]
    res = lsim(Gcl, disturbance, 100)
    f3 = plot(res, ylabel=["y" "u"]); ylims!((-0.05, 0.3), sp = 1)
    f4 = gangoffourplot(prob, w, lab = "", legend = :bottomright)

    plot(f1, f2, f3, f4, titlefontsize=10)
end

plot(plots..., size=(1200,1000))



# ==============================================================================
## LinearMPC
# The example below illustrate how we can convert the LQGProblem into an MPC contorller by loading LinearMPC.jl
# We then perform a rather low-level simulation with a manual loop, where we form the observer `obs` and step the plant
# ==============================================================================

using LinearMPC

Gd = add_low_frequency_disturbance(G, ϵ = 1e-6, measurement=false)

C  = Gd.C
Q1 = 100diagm([1.0])    # output cost matrix
R1 = diagm([0.001, 1])  # Dynamics noise covariance
R2 = I(ny)              # Measurement noise covariance
Gde = ExtendedStateSpace(Gd, B1=I) # Since B1=I, R1 has size determined by state dimension, but C1=C, so Q1 has size determined by the output dimension
prob = LQGProblem(Gde, Q1, Q2, R1, R2)
obs = observer_predictor(prob, direct=false) # If a predictor is used, the observer update should be carried out in the end of the loop, if we use the filter below, we should instead perform the observer update in the beginning of the loop directly after obtaining the new measurement but before computing a new control signal.
obs = observer_filter(prob)

mpc = LinearMPC.MPC(prob; N=20, umin=[-3], umax=[3])
x = zeros(G.nx)     # True plant state
xh = zeros(Gd.nx)   # Observer state
X = [x[]]           # For storage
U = Float64[]
u_mpc = [0.0]       
for i = 1:50
    y = G.C*x                          # Compute the true measurement output
    xh = obs.A * xh + obs.B*[u_mpc; y] # Predict one step with the observer, u here is the control signal from the previous iteration, if using the predictor, use u from the current iteration and perform the observer update in the end of the loop instead
    u_disturbance = i * Ts ≥ 10 ? 1.0 : 0.0
    r = [1.0] # output reference
    u_mpc = LinearMPC.compute_control(mpc, xh; r) # Call MPC optimizer with estimated state
    u_tot = u_mpc .+ u_disturbance # Total input is control signal + disturbance
    x = G.A*x + G.B*u_tot   # Advance the true plant state
    push!(X, x[])           # Store data for plotting
    push!(U, u_mpc[])
end


plot(X*G.C', layout=2, sp=1, label="\$y\$")
plot!(U, sp=2, label="\$u\$")
hline!([1.0 3.0], linestyle=:dash, color=:black, label=["Reference" "\$u_{max}\$"], sp=[1 2])


using Test
@test (G.C*X[end])[] ≈ 1 rtol=1e-4