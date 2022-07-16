# Mixed-sensitivity ``H_\infty`` control design

In this tutorial, we will design a controller for a DC-servo (electrical motor).

The servo takes a torque command as input and we are interested in controlling the angle of the output shaft. An approximate model for a DC servo, valid for low frequencies, is

$\tau = k ( J \ddot{\phi} + c\dot{\phi} )$

where $\tau$ is the torque and $\phi$ is the angle. The parameter $J$ denotes the inertia of the motor and load, $c$ is a viscous friction parameter and $k$ is a torque constant (gain). 

This is a simple SISO example with a pole in the origin. Poles on the stability boundary are problematic for several numerical routines, and this tutorial demonstrates a common trick applied in these situations. 

We start by defining the process model. We will use the parameter values $J=1, c=0.12, k=11.2$ which yields the transfer function

$G(s) = \dfrac{k}{s(Js + c)} = \dfrac{11.2}{s^2 + 0.12s + 0}$

```@example hinfdesign
using RobustAndOptimalControl, ControlSystems, Plots
Gtrue = tf([11.2], [1, 0.12, 0])
nothing # hide
```

When designing a controller using $H_\infty$ synthesis, we formally specify an optimization problem where the cost function is the $H_\infty$ norm of a suitably chosen transfer function, and the optimization variable is the controller. The function [`hinfpartition`](@ref) helps us to build the transfer function that appears in the following $H_\infty$ optimization problem

```math
\operatorname{minimize}_K \begin{Vmatrix}
W_S S \\
W_U KS \\
W_T T
\end{Vmatrix}_\infty
```


where $K$ is the controller and $S$ is the sensitivity function $(1+GK)^{-1}$. The transfer functions $W_S$, $W_T$ and $W_U$ are weight functions that emphasize different frequency ranges, for example, if $W_S(i\omega)$ is large for a particular frequency $\omega$, then $S$ is forced to be small at $\omega$ in order for the $H_\infty$ norm to be small. 

In this case, we will omit the penalty on the transfer function ``T``, and thus solve
```math
\operatorname{minimize}_K \begin{Vmatrix}
W_S S \\
W_U KS
\end{Vmatrix}_\infty
```

```@example hinfdesign
# Sensitivity weight function
WS = makeweight(1e5, 0.05, 0.5) |> tf

# Output sensitivity weight function. Increase this value to penalize controller effort more
WU = ss(1.0)

# Complementary sensitivity weight function
WT = [] # We do not put any weight on T in this example
nothing # hide
```

To solve this problem, we partition the system $P$ in a two-input, two output configuration such that 
$\operatorname{lft}_{l}(P,K)$ forms the system we want to minimize the $H_\infty$ norm of with respect to ``K``. To help with this partitioning, we have the function [`hinfpartition`](@ref):
```@example hinfdesign
P = hinfpartition(Gtrue, WS, WU, WT)
nothing # hide
```

The object `P` will now be of type [`ExtendedStateSpace`](@ref), and represent the following `P`, with two input ports `w,u` and two output ports `z,y`:
```
     ┌─────────┐
z◄───┤         │◄────w
     │    P    │
y┌───┤         │◄───┐u
 │   └─────────┘    │
 │                  │
 │      ┌───┐       │
 │      │   │       │
 └─────►│ K ├───────┘
        │   │
        └───┘
```
The operation `lft(P, K)` forms the feedback interconnection in the diagram, and the $H_\infty$ optimization will minimize the $H_\infty$ norm of the transfer function from ``w`` to ``z`` with respect to ``K``.

Before solving, we may check if the synthesis problem is feasible
```@example hinfdesign
hinfassumptions(P, verbose=true)
```
The problem is not feasible, in this case due to the integrator pole on the stability margin. We thus modify the plant description by moving the integrator pole in the origin slightly ($\vareplsilon$) into the stable region

```@example hinfdesign
ε = 1e-4
G = tf([11.2], [1, 0.12]) * tf([1], [1, ε])
nothing # hide
```
If we now perform the assumption check again, it passes
```@example hinfdesign
P = hinfpartition(G, WS, WU, WT)
hinfassumptions(P)
```

## Synthesize the H-infinity optimal controller
With the problem properly defined, we may call [`hinfsynthesize`](@ref) to solve it. The result is the controller ``K`` and a performance index $\gamma$:
```@example hinfdesign
K, γ = hinfsynthesize(P, γrel=1.05)
γ
```
The achieved performance level is indicated by $\gamma$, this number is the norm we are optimizing (the $H_\infty$ norm) and it should be as low as possible.

!!! note "$H_2$ optimization"
    If we instead of [`hinfsynthesize`](@ref) had called [`h2synthesize`](@ref), we had solved the same optimization problem, but under the ``H_2`` norm instead of the ``H_\infty`` norm.

For verification purposes, we may extract some transfer functions defining common sensitivity functions:
```@example hinfdesign
Pcl, S, KS, T = hinfsignals(P, G, K)
nothing # hide
```
We may also verify that the closed-loop system has $H_\infty$ norm $\gamma$
```@example hinfdesign
Pcl == lft(P, K)
isapprox(hinfnorm2(Pcl)[1], γ, rtol=1e-5)
```

## Plot the specifications
The resulting sensitivity functions can be plotted together with the inverse weighting functions multiplied by $\gamma$ to get a feeling for where the constraints are active. 
```@example hinfdesign
specificationplot([S, KS, T], [WS, WU, WT], γ)
```
In this case, the noise amplification constraint is active between about $10^0 - 10^{1}$ rad/s and the sensitivity function is pushed down for lower frequencies. In this case, the complimentary sensitivity function $T = (1+GK)^{-1}GK$ has desireable properties without us specifying and penalizing this function in the optimization problem. If we would like to push this function down at certain frequencies, we may specify a non-empty weight $W_T$ as well. The transfer function ``KS`` (labeled `CS` in the figure) has some natural roll-off for high frequencies. If we want steeper roll-off (more filtering), we could change the weight function `WU` to a transfer function with high gain for high frequencies.


## Generate C-code for the controller
The controller `K` is given in the form of a statespace system. The package [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl) can be used to generate C-code for such systems, making it easy to implement the controller.