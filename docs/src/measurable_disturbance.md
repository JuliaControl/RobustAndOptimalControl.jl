# Feedforward from known disturbances

This example will demonstrate how you can make use of a known/measurable disturbance and achieve effective disturbance rejection using an $\mathcal{H}_2$ controller. For simplicity, we will consider a simple first-order system $G$

```math
\begin{aligned}
\dot{x} &= -ax + b(u + d) \\
y_x &= cx + e_y \\
y_d &= d + e_d
\end{aligned}
```
where a load disturbance $d$ is acting on the input of the system. We have access to a measurement ``y_x`` of the state and a measurement ``y_d`` of the disturbance. Both of them are assumed corrupted by measurement noise (this fact is crucial, if there is no measurement noise on the disturbance measurement, the controller will ignore the measurement of the state).

We start by defining the process model.

```@example LQG_MEASURABLE_DIST
using RobustAndOptimalControl, ControlSystems, Plots, LinearAlgebra
G = ss(tf(1, [10, 1])) # Process model
```

We then choose the parameters of the $\mathcal{H}_2$ controller, i.e., noise variances and how much we penalize state and inputs. When tuning a $\mathcal{H}_2$ controller, we do not choose matrices in quite the same way as when we tune an LQG controller, but an LQG controller is actually a $\mathcal{H}_2$ in disguise, see sec. 9.3.3 in Skogestad for more details. 

We will partition the system according to
```math
\begin{bmatrix}
\dot x \\ \hline z_y \\ z_u \\ \hline y_x \\ y_d
\end{bmatrix}=
\begin{bmatrix}
A   & B_d & 0 & 0 & B_u \\ \hline
C_z & 0   & 0 & 0 & 0   \\
0   & 0   & 0 & 0 & r   \\ \hline
C   & 0   & 1 & 0 & 0   \\
0   & 1   & 0 & 1 & 0   \\
\end{bmatrix}
\begin{bmatrix}
x \\ \hline d \\ e_y \\ e_d \\ \hline u
\end{bmatrix}
```


where ``C_z = qC`` and ``B_d = w_d B``.

```@example LQG_MEASURABLE_DIST
σye = 1  # standard deviation of output measurement
σyd = 10 # standard deviation of known-disturbance measurement
wd = 100 # Disturbance suppression weight

q = 10   # Penalty on output
r = 0.1  # Penalty on input

Ge = ExtendedStateSpace(G,
    B1 = [wd*G.B 0 0],        # Load disturbance, y measurement noise, d measurement noise
    C1 = [q*G.C; 0G.C],    # Q1 = C1'C1
    C2 = [G.C; 0G.C],      # Measure output and load disturbance, load disturbance is not a function of state
    D12 = [0; r],              # Penalize control action Q2 = D12'D12
    D21 = [0 σye 0; 1 0 σyd],  # direct feedthrough of load disturbance and measurement noise
)

h2 = false
if h2
    K, Cl = h2synthesize(Ge, 100)
else
    K, γ, mats = hinfsynthesize(Ge)
    Cl = lft(Ge, K)
end
```

```@example LQG_MEASURABLE_DIST
x0  = zeros(G.nx) # Initial condition
disturbance = (x, t) -> [10 < t * Ts < 20; 0randn(); 0randn()] # This is our load disturbance, a step at ``t = 10``

res = lsim(Cl, disturbance, 100, dt=0.001, adaptive=false)
plot(res)
```

Before we feel confident about deploying the LQG controller, we investigate its closed-loop properties. The (negative) feedback controller is given by `-K[1,1]` (the - sign since [`h2synthesize`](@ref) returns a positive-feedback controller).

```@example LQG_MEASURABLE_DIST
w = exp10.(LinRange(-3, 4, 300))
gangoffourplot(G, -K[1,1], w, lab = "", legend = :bottomright)
```

We see that our design led to a system with reasonable disturbance-rejection properties, however, since we assume that we can measure the disturbance, it's relevant to also consider the transfer function from this disturbance to the output.

This transfer function is given by `Cl[1,1]` and has been scaled by the disturbance-suppression weight ``w_d`` and the performance penalty weight ``q``, so we rescale it to get back the original units.
```@example LQG_MEASURABLE_DIST
bodeplot!((1/q/wd)Cl[1,1], w, plotphase=false, sp=2, lab="Gyd", l=:dash)
```
we see that it overlaps completely with the classical sensitivity function `PS`, which is expected.