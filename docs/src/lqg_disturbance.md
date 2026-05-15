# Disturbance modeling and rejection with LQG controllers

This example will demonstrate how you can add disturbance models to ta plant model and achieve effective disturbance rejection using an LQG controller. For simplicity, we will consider a simple first-order system $G$

```math
\begin{aligned}
\dot{x} &= -ax + b(u + d) \\
y &= cx
\end{aligned}
```

where a load disturbance $d$ is acting on the input of the system. This is a simple and very common model for load disturbances. In this example, we will let $d$ be a unit step at time $t=10$, this will effectively create a LQG controller with *integral action*.

We will begin by setting up the LQG problem and solve it without andy disturbance model. For details regarding the setup of an LQG problem, see, the [`LQGProblem`](@ref) documentation.

We start by defining the process model and discretize it using zero-order hold.

```@example LQG_DIST
using RobustAndOptimalControl, ControlSystemsBase, Plots, LinearAlgebra
Ts = 1 # Sample time
G = c2d(ss(tf(1, [10, 1])), Ts) # Process model
```

We then choose the parameters of the LQG controller, i.e., the cost matrices ``Q_1, Q_2`` as well as the covariance matrices ``R_1, R_2``
```@example LQG_DIST
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
```

As we can see, the controller appears to do next to nothing to suppress the disturbance. The problem is that the Kalman filter does not have a model for such a disturbance, and its estimate of the state will thus be severely biased.

The next step is to add a disturbance model to the plant model. Since the disturbance if of low-frequency character (indeed, its transfer function is $1/s$), we make use of the function `add_low_frequency_disturbance`

```@example LQG_DIST
Gd = add_low_frequency_disturbance(G, ϵ = 1e-6) # The ϵ moves the integrator pole slightly into the stable region
nx = Gd.nx
```

There is no point trying to penalize the disturbance state in the LQR problem, it's not controllable, we thus penalize the output only, which we can write as

```math
y^T Q_1 y = (Cx)^T Q_1 Cx = x^T (C^T Q_1C) x
```

```@example LQG_DIST
C  = Gd.C
Q1 = 100C'diagm(ones(G.nx)) * C # state cost matrix
x0 = zeros(nx)
nothing # hide
```

We also provide new covariance matrices for the Kalman filter where the entry of the state-covariance matrix that corresponds to the disturbance state (the second and last state) determines how fast the Kalman filter integrates the disturbance. We choose a large value (1), implying fast integration

```@example LQG_DIST
R1 = diagm([0.001, 1])
R2 = I(ny)
SQ = [zeros(nx-nu, nu); Q2] # Adding Q2 to the x-u cross term achieves perfect integral action without manually modifying the computed feedback gain, if not all inputs are augmented with integral action, we should add Q2[augmented_inds, :] instead
prob = LQGProblem(Gd, Q1, Q2, R1, R2; SQ=SQ)
Gcl  = [G_PS(prob); -comp_sensitivity(prob)] # -comp_sensitivity(prob) is the same as the transfer function from load disturbance to control signal
res  = lsim(Gcl, disturbance, 100)
plot(res, ylabel=["y" "u"]); ylims!((-0.05, 0.3), sp = 1)
```

This time, we see that the controller indeed rejects the disturbance and the control signal settles on -1 which is exactly what's required to counteract the load disturbance of +1.

```@example LQG_DIST
using Test
@test lqr(prob)[end] ≈ 1.0 # The gain from estimated disturbance state to control signal should be 1.0, indicating perfect integral action
```

Before we feel confident about deploying the LQG controller, we investigate its closed-loop properties.

```@example LQG_DIST
w = exp10.(LinRange(-3, log10(pi / Ts), 200))
gangoffourplot(prob, w, lab = "", legend = :bottomright)
```

We see that our design led to a system with a rather high peak in sensitivity. This is an indication that we perhaps added too much "integral action" by a too fast observer pole related to the disturbance state. Let's see how a slightly more conservative design fares:

```@example LQG_DIST
R1 = diagm([0.001, 0.2]) # Reduce the noise on the integrator state from 1 to 0.2
R2 = I(ny)
prob = LQGProblem(Gd, Q1, Q2, R1, R2)

Gcl = [G_PS(prob); -comp_sensitivity(prob)]
res = lsim(Gcl, disturbance, 100)
f1 = plot(res, ylabel=["y" "u"]); ylims!((-0.05, 0.3), sp = 1)
f2 = gangoffourplot(prob, w, lab = "", legend = :bottomright)

plot(f1, f2, titlefontsize=10)
```

We see that we now have a slightly larger disturbance response than before, but in exchange, we lowered the peak sensitivity and complimentary sensitivity from (1.51, 1.25) to (1.31, 1.11), a more robust design. We also reduced the amplification of measurement noise ($CS = C/(1+PC)$). To be really happy with the design, we should probably add high-frequency roll-off as well.

## Extracting the controller

Above we simulated the closed-loop response by building closed-loop transfer functions directly with [`G_PS`](@ref) and [`comp_sensitivity`](@ref). To obtain the controller used under the hood, call [`observer_controller`](@ref) for the pure feedback controller, or [`extended_controller`](@ref) for a 2-DOF controller that takes a separate state reference.

## Explicit error integration via LQI

The disturbance-model approach above achieves integral action implicitly: the Kalman filter estimates a slow disturbance state, and the LQR gain on that state acts as an integrator. A more direct alternative is to design a Linear-Quadratic-Integral (LQI) controller, in which output-error integrators are augmented into the plant state and the LQR cost penalizes the integrated error explicitly. The interface is [`lqi`](@ref) (gain only) and [`lqi_controller`](@ref) (full controller that takes references and measurements as inputs).

We re-use the original (un-augmented) plant `G` from the top of this page:

```@example LQG_DIST
nx = G.nx
nu = G.nu
ny = G.ny

# Cost on the augmented state x_a = [x; ∫e].
# The trailing diagonal block weights the integral of the tracking error;
# increasing it makes the controller act more aggressively on accumulated error.
Q1 = cat(100.0*I(nx), 100*I(ny); dims=(1, 2))
Q2 = 0.01*I(nu)

# Observer for the un-augmented plant
R1 = 0.001I(nx)
R2 = I(ny)
K   = kalman(G, R1, R2)
obs = observer_predictor(G, K; output_state = true)

C = lqi_controller(G, obs, Q1, Q2)   # controller with inputs [r; y] and output u
```

To inject the load disturbance from earlier (a unit step added at the plant input), we close the loop using only the measurement column of `C` and feed the disturbance into the plant input:

```@example LQG_DIST
Cy    = C[:, :y_plant]         # feedback channel (r=0)
Kctrl = -Cy                    # equivalent negative-feedback controller
Gcl_y = feedback(G, Kctrl)     # disturbance → y
Gcl_u = -Kctrl * Gcl_y         # disturbance → u
Gcl   = [Gcl_y; Gcl_u]
res = lsim(Gcl, disturbance, 100)
@test res.y[:, end] ≈ [0, -1] atol=1e-3
plot(res, ylabel = ["y" "u"]); ylims!((-0.05, 0.3), sp = 1)
```

The control signal again settles on `-1`, exactly counteracting the load disturbance. Compared to the disturbance-model design, the LQI formulation lets us tune the strength of the integral action directly through the augmented-state cost block in `Q1`, rather than indirectly through the disturbance-state noise covariance in `R1`.
