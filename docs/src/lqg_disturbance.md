# Disturbance modeling and rejection with LQG controllers

A constant unmeasured input disturbance can be estimated by augmenting the observer
model with constant disturbance variables. The physical feedback gain is designed
on the plant, and the disturbance estimate is subtracted from the control input.

For the scalar plant
```math
\dot{x} = -ax + b(u+d), \qquad y = cx,
```
the augmented observer includes a variable with dynamics ``\dot{d}=0`` in continuous
time, or ``d^+=d`` in discrete time. This variable is uncontrollable and is excluded
from the state-feedback Riccati equation.

## Physical plant and baseline design

```@example LQG_DIST
using RobustAndOptimalControl, ControlSystemsBase, Plots, LinearAlgebra, Test
Ts = 1
G = c2d(ss(tf(1, [10, 1])), Ts)
nx, nu, ny = G.nx, G.nu, G.ny
Q1 = 100I(ny)
Q2 = 0.01I(nu)
R1 = 0.001I(nx)
R2 = I(ny)
Pe = ExtendedStateSpace(G; B1=I, C1=G.C)
prob = LQGProblem(Pe, Q1, Q2, R1, R2)
disturbance = (x, t) -> t * Ts ≥ 10
plot(lsim(G_PS(prob), disturbance, 100))
```

The baseline observer does not model the constant disturbance, so its estimate is
biased under a sustained load.

## Exact constant-disturbance model

```@example LQG_DIST
Gd = add_low_frequency_disturbance(G)
Pe = ExtendedStateSpace(Gd; B1=I, C1=Gd.C)
disturbance_indices = collect((G.nx+1):Gd.nx)
disturbance_input_indices = collect(1:nu)
R1 = Diagonal([0.001, 1.0])
prob = LQGProblem(Pe, Q1, Q2, R1, R2;
    disturbance_indices, disturbance_input_indices)

L = lqr(prob)
@test L[:, disturbance_indices] ≈ I(nu)
@test isstable(feedback(G, observer_controller(prob)))
@test norm(dcgain(G_PS(prob))) < 1e-8
Gcl = [G_PS(prob); -input_comp_sensitivity(prob)]
res = lsim(Gcl, disturbance, 100)
plot(res, ylabel=["y" "u"])
```

`disturbance_indices` refers to augmented state-variable indices, whereas
`disturbance_input_indices` refers to control-input indices. Their ordering defines
the pairing. For a subset of inputs, pass their indices and augment the observer
using the corresponding columns of `G.B`. Model reduction must precede this
augmentation so that the disturbance coordinates remain identifiable.

If the input indices are omitted, the input mapping is inferred when the physical
input matrix has full column rank and exactly represents the disturbance coupling.
Redundant actuator configurations require an explicit input pairing. Unmatched
disturbances and direct disturbance contributions to outputs are unsupported by
this interface.

The Kalman filter uses the complete augmented model. The covariance of each
disturbance variable determines the estimator response to changing disturbances.
The controller implements ``u=-L_x\hat{x}-\hat{d}``; its disturbance compensation
gain is specified independently of the physical LQR optimization.

```@example LQG_DIST
w = exp10.(LinRange(-3, log10(pi / Ts), 200))
gangoffourplot(prob, w, lab="", legend=:bottomright)
```

A smaller disturbance covariance reduces the estimator response rate and changes
the trade-off between disturbance rejection and measurement-noise amplification.

```@example LQG_DIST
R1 = Diagonal([0.001, 0.2])
prob = LQGProblem(Pe, Q1, Q2, R1, R2;
    disturbance_indices, disturbance_input_indices)
Gcl = [G_PS(prob); -input_comp_sensitivity(prob)]
f1 = plot(lsim(Gcl, disturbance, 100), ylabel=["y" "u"])
f2 = gangoffourplot(prob, w, lab="", legend=:bottomright)
plot(f1, f2)
```

## Controller and reference response

`observer_controller(prob)` returns the feedback controller with the standard
negative-feedback convention. `extended_controller(prob)` provides state-reference
and measurement inputs with the feedback sign included in the controller.

The physical plant is returned by `system_mapping(prob)`. The augmented model is
available through `prob.sys`; `performance_mapping(prob)` and `lft(prob)` retain
its disturbance channels. Physical closed-loop stability must be assessed without
including the autonomous disturbance generator as a plant mode.

```@example LQG_DIST
Cfb = observer_controller(prob)
Cff = ff_controller(prob)
Gref = closedloop(prob) * static_gain_compensation(prob)
@test isstable(Gref)
@test dcgain(Gref) ≈ I(ny)
plot(step(Gref, 100))
```

`G_PS(prob)*Cff` has the same reference transfer function but may contain canceling
controller modes in its realization. `Gref` avoids these redundant modes.

## Explicit error integration via LQI

The disturbance-model approach above achieves integral action implicitly: the Kalman filter estimates a constant disturbance, and the controller compensates for that estimate. A more direct alternative is to design a Linear-Quadratic-Integral (LQI) controller, in which output-error integrators are augmented into the plant state and the LQR cost penalizes the integrated error explicitly. The interface is [`lqi`](@ref) (gain only) and [`lqi_controller`](@ref) (full controller that takes references and measurements as inputs).

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

To inject the load disturbance from earlier (a unit step added at the plant input), we close the loop with the advanced [`feedback`](@ref) interface on named systems. Two features of that interface keep the wiring direct:
- listing `:u_plant` in both `w1` and `u1` makes the plant's input serve simultaneously as the feedback control input from `C` and as the external summing point for the load disturbance, so no extra disturbance port has to be added to `G`;
- omitting the reference `:y_plant_r` from both `w2` and `u2` is equivalent to `r = 0` without exposing it as a closed-loop input, so the original scalar `disturbance` function can be passed to `lsim`.

```@example LQG_DIST
G_named = named_ss(G, u = :u_plant, y = :y_plant)
Gcl = feedback(G_named, C;
        w1 = :u_plant,            # external input: load disturbance at the plant input
        u1 = :u_plant,            # feedback path drives the same physical input → both sum
        z1 = :y_plant,            # keep plant output
        # w2 = :y_plant_r,        # uncomment to also expose the reference as a closed-loop input
        u2 = :y_plant,            # C feedback input ← plant output
        z2 = :y_L,                # keep controller output (the control signal u)
        pos_feedback = true)      # `lqi_controller` already bakes in the negative-feedback sign
res = lsim(Gcl, disturbance, 100)
@test res.y[:, end] ≈ [0, -1] atol=1e-3
plot(res, ylabel = ["y" "u"]); ylims!((-0.05, 0.3), sp = 1)
```

The control signal again settles on `-1`, exactly counteracting the load disturbance. Compared to the disturbance-model design, the LQI formulation lets us tune the strength of the integral action directly through the augmented-state cost block in `Q1`, rather than indirectly through the disturbance-state noise covariance in `R1`.


## 2-DOF tracking with [`extended_controller`](@ref)

The reference-signal response of an LQG controller can be improved without adding integral action by proper reference feedforward. [`extended_controller`](@ref) returns an [`ExtendedStateSpace`](@ref) controller with inputs `[xᵣ; y]` and output `u`. The reference enters the observer dynamics through `B1 = (B − KD)·L`, so the observer estimate `x̂` tracks the true state `x` even when `xᵣ ≠ 0` (in contrast to a pure feedforward path that bypasses the observer).

We re-use the original plant `G` and build a fresh LQG problem (no disturbance model needed here):

```@example LQG_DIST
Q1 = 100I(G.nx)
Q2 = 0.01I(G.nu)
R1 = 0.001I(G.nx)
R2 = I(G.ny)
prob_2dof = LQGProblem(G, Q1, Q2, R1, R2)

# z=[1] returns the closed-loop transfer function from xᵣ to plant output 1 as a
# second return value, useful for computing DC-gain compensation.
Ce, cl_xr_to_y = extended_controller(prob_2dof, z=[1])
```

The closed-loop DC gain from state reference to output is not unity in general — for `xᵣ = x_ss` to hold, the plant would need to contain an integrator. For this stable first-order plant we use the second return value to compute a static pre-compensation:

```@example LQG_DIST
gain_comp = inv(dcgain(cl_xr_to_y)[1])
res = lsim(gain_comp * cl_xr_to_y, (x,t) -> min(t/10, 1), 20)
@test res.y[end] ≈ 1 atol=1e-2
plot(res, ylabel="y")
```

The step response settles on `1`, confirming that the pre-compensated 2-DOF controller tracks unit references at DC. For plants with an integrator (e.g., a cart-position channel from a velocity-controlled actuator) `dcgain(cl_xr_to_y)` is already `1` and no compensation is needed.

Note, without the integral action in the controller, any error in the DC gain of the model would still lead to a steady-state error in the reference-signal response.


## Integral action removes the steady-state error under model mismatch

To make the point of the previous paragraph concrete, we apply both controllers — the integrating `C` from [`lqi_controller`](@ref) and the non-integrating `Ce` from [`extended_controller`](@ref), each designed earlier on the nominal plant `G` — to a *perturbed* plant whose gain is 20 % larger than the design model:

```@example LQG_DIST
G_pert   = c2d(ss(tf(1.2, [10, 1])), Ts)   # 20 % larger gain than the design model G
G_pert_n = named_ss(G_pert)

# Close the reference-tracking loop directly on the named systems
#   w1 = :y_plant_r  → the reference is the external input
#   u1 = :y_plant    → wire C's feedback input to the plant output
#   z2 = G_pert_n.y  → keep only the plant output (suppress C's u as an output)
#   pos_feedback = true → lqi_controller already bakes in the negative sign of the feedback path
Gcl_pert_int = feedback(C, G_pert_n;
                    w1 = :y_plant_r,
                    z1 = Symbol[],
                    z2 = G_pert_n.y,
                    u1 = :y_plant,
                    pos_feedback = true)

# The controller includes the negative measurement-feedback sign.
Ce_named = named_ss(ss(Ce), y = C.y, u = C.u)
Gcl_pert_noint = feedback(G_pert_n, Ce_named;
                    w1 = Symbol[],
                    w2 = :y_plant_r,
                    u2 = :y_plant,
                    z1 = G_pert_n.y,
                    z2 = Symbol[],
                    pos_feedback = true)

res_int   = lsim(Gcl_pert_int,               (x,t)->min(t/10, 1), 60)
res_noint = lsim(gain_comp * Gcl_pert_noint, (x,t)->min(t/10, 1), 60)
@test res_int.y[end]   ≈ 1   atol=1e-3
@test res_noint.y[end] ≈ dcgain(gain_comp * Gcl_pert_noint)[] atol=1e-2
@test !isapprox(res_noint.y[end], 1; atol=1e-3)
plot(res_int, lab = "lqi_controller")
plot!(res_noint, ylabel = "y", lab = "extended_controller")
```

The integrating controller rejects the steady-state error caused by this gain mismatch. The non-integrating controller has a nonzero reference error whose magnitude depends on the plant mismatch and the feedback design.
