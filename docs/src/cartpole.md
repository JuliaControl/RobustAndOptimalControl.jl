# Control design for a pendulum on a cart
In this example we will consider control design for the basic inverted pendulum on a cart. This system has two equilibria, one where the pendulum is hanging straight down, and one where it's balancing straight up. The upper one is unstable, making it slightly more interesting to design a controller for (even if the lower equilibrium is highly relevant, it's a good model for an overhead crane moving goods).

## System model
In this tutorial, we assume that we have the nonlinear dynamics of the system encodeed as a julia function `ẋ = cartpole(x, u)`, and linearize this to get a statespace system
```math
\begin{aligned}
ẋ &= Ax + Bu\\
y &= Cx
\end{aligned}
```
We make use of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl/) for the linearization. We start by defining the dynamics function
```@example PENDCART
using ControlSystems, RobustAndOptimalControl, ForwardDiff, LinearAlgebra, Plots
default(label="") # hide

function cartpole(x, u)
    mc, mp, l, g = 1.0, 0.2, 0.5, 9.81

    q  = x[1:2]
    qd = x[3:4]

    s = sin(q[2])
    c = cos(q[2])

    H = [mc+mp mp*l*c; mp*l*c mp*l^2]
    C = [0.1 -mp*qd[2]*l*s; 0 0]
    G = [0, mp * g * l * s]
    B = [1, 0]

    qdd = -H \ (C * qd + G - B * u[1])
    return [qd; qdd]
end

nu = 1    # number of control inputs
nx = 4    # number of states
ny = 2    # number of outputs (here we assume that the cart position and the pendulum angle are measurable)
nothing # hide
```

## Linearization
The next step is to choose an operating point around which to linearize and to calculate the Jacobians $A$ and $B$:
```@example PENDCART
x0 = [0, π, 0, 0]
u0 = [0]

Ac = ForwardDiff.jacobian(x->cartpole(x, u0), x0)
Bc = ForwardDiff.jacobian(u->cartpole(x0, u), u0)
Cc = [1 0 0 0; 0 1 0 0]
Λ = Diagonal([0.4, deg2rad(25)]) # Maximum output ranges
Cc = Λ\Cc # This normalizes expected outputs to be ∈ [-1, 1], a good practice for MIMO systems
nothing # hide
```
we package everything into a [`StateSpace`](@ref) object and visualize its poles and zeros:
```@example PENDCART
sys = ss(Ac, Bc, Cc, 0)
```

```@example PENDCART
pzmap(sys)
```

## Control design

We will design a number of different controllers. We will start with a basic PID controller. Since the PID controller in its standard form really only handles SISO systems, we will also design a state-feedback controller with an observer to estimate the full state vector $x$ based on the two measurements $y$. Lastly, we will attempt to "robustify" the state-feedback controller using the [`glover_mcfarlane`](@ref) procedure.

Since the system has an unstable pole $p \approx 4.85$rad/s, there wil be fundamental limitations on the performance of the closed loop system. A common rule-of-thumb (see, e.g., Åström and Murray) is that a single RHP pole $p$ puts a *lower* limit on the gain crossover frequency $\omega_{gc} > 2p$, something to take into consideration when tuning our controllers. 

## PID controller
Since the PID controller only accepts a single measurement, we choose the measurement of the pendulum angle for feedback. While doing so, we notice that the number of states in the model can be reduced by the function [`sminreal`](https://juliacontrol.github.io/ControlSystems.jl/latest/lib/synthesis/#ControlSystems.sminreal-Tuple{StateSpace})
```@example PENDCART
P = sminreal(sys[2,1]) # Position state goes away, not observable
```
this indicates that the state corresponding to the position of the cart is not observable from the measurement of the pendulum angle. This is slightly worrisome, but we nevertheless proceed to design a controller.
By using a single measurement only, we have also introduced a zero in the system
```@example PENDCART
pzmap(P)
```
A PID controller can be constructed using the function [`pid`](https://juliacontrol.github.io/ControlSystems.jl/latest/lib/synthesis/#ControlSystems.pid-Tuple{}). We start our tuning by a simple P controller
```@example PENDCART
C = pid(1, 0, 0)
```
We will attempt to perform loop shaping using the PID controller, and plot the stability margins in a Bode plot using the function `marginplot`
```@example PENDCART
w = exp10.(LinRange(-2.5, 3, 500))
function pid_marginplot(C)
    f1 = marginplot(P*C, w)
    vline!([2*4.85], sp=1, lab="Fundamental limitation", l=(:dash, :black))
    ylims!((1e-3, 1e2), sp=1)
    f2 = nyquistplot(P*C)
    plot(f1, f2)
end
pid_marginplot(C)
```
We notice that the gain of the loop-transfer function $L = PC$ is much too low, and increase it, we also notice that the Nyquist plot fails to encircle to critical point, which it has to do once since we have one unstable pole. We will solve this in the end by adding integral action, but proceed for now to shape other parts of the loop. We start by lifting the Bode curve by increasing the gain:
```@example PENDCART
C = pid(20, 0, 0)
pid_marginplot(C)
```
we are now getting close to the rule-of-thumb for $\omega_{gc}$, but have a low loop gain at low frequencies. Remember, to get good disturbance rejection, we typically want a high loop gain at low frequencies. We also have an extremely small phase margin at 0.66 degrees. To fix the phase margin, we add some derivative gain. While adding derivative gain, it's also a good idea to add noise filtering (with a pure derivative term, the PID controller is not proper and can not be realized as a statespace system)
```@example PENDCART
C = pid(20, 0, 0.2, Tf=0.01)
pid_marginplot(C)
```
The derivative term lifted the phase at $\omega_{gc}$ and we now have very nice phase margins. We also got a slight increase in $\omega_{gc}$ while at it. 

The closed-loop system will still be unstable since the Nyquist curve fails to encircle the point -1, something we can check by calling
```@example PENDCART
isstable(feedback(P*C))
```
We make the Nyquist curve wrap around the -1 point by adding integral gain:
```@example PENDCART
C = pid(20, 1.25, 0.2, Tf=0.01)
pid_marginplot(C)
```
Now, the Nyquist curve looks fine and the system is stable
```@example PENDCART
isstable(minreal(feedback(P*C)))
```


If we simulate a disturbance acting on this system (`feedback(P, C)` is the transfer function from load disturbance to output)
```@example PENDCART
plot(step(feedback(P,C), 8), ylab="ϕ")
```
we see that we have a reasonable disturbance response. 

To verify robustness properties, we plot the gang-of-four sensitivity functions:
```@example PENDCART
f1 = gangoffourplot(P,C,w, Ms_lines=[1.4], Mt_lines=[1.5])
f2 = nyquistplot(P*C, Ms_circles=[1.4], Mt_circles=[1.5], ylims=(-2, 2), xlims=(-4,1))
plot(f1, f2, size=(1000,800))
```
This all looks nice and we appear to have reasonable robustness margins, the Nyquist curve stays outside the $M_S = 1.4$ circle and the $M_T = 1.5$ circle.

However, there is a dragon lurking behind these plots. Remember the state corresponding the the cart position that was removed above? What has happened to this state? To investigate this, we form an [`ExtendedStateSpace`](@ref) model where we have both cart position and pendulum angle as controlled outputs, while keeping only the pendulum angle as measured output:
```@example PENDCART
Pe = ExtendedStateSpace(sys, C2 = sys.C[2:2, :]) # Indicate that we can only measure the pendulum angle
Gecl = feedback(Pe, ss(C)) |> minreal
plot(step(Gecl, 8), ylab=["Cart pos" "ϕ"])
```



We see that the cart position drifts away without ever thinking about stopping. Indeed, the PID controller is unaware of this and can not really do anything about it. We could attempt to design a second control loop that would close the loop around the cart position, but we would have to carefully manage the interactions between the two loops. Instead, we move on to a state-feedback design, a methodology that makes handling multiple outputs much more straightforward. 

## Pole placement and observer design
The design of a state-feedback controller typically involves two steps, designing the feedback gain and designing an observer. We will arrive at the feedback gain through pole placement, but will design the observer as a Kalman filter, i.e., by solving a Riccati equation rather than using Ackermann's formula. 

When performing pole placement, there are a number of design guidlines that help you arrive at a robust design. One of these are that past process poles should be matched with an equally fast closed-loop pole. We can get an overview of the open-loop poles with `dampreport`
```@example PENDCART
dampreport(sys)
```
we see that we have two poles at roughly $\pm 4.85$rad/s, and almost two integrators. We thus keep the fast pole, and place the unstable pole at the same location (same bandwidth but stable instead of unstable). We also try to move the integrator poles to -5 to make the system nice and fast. 

```@example PENDCART
desired_poles = [-4.85, -4.85, -5, -5]
L = place(sys, desired_poles, :c)
```

For the observer, we make use of the function `kalman`. We choose the covariance matrices `R1, R2` that determine the amount of noise acting on the system and on the measurements respectively. We assume that there are two noise components, both entering as forces. One disturbance force acts on the cart and the other on the pendulum. We indicate this using the matrix ``B_w``. 
```@example PENDCART
Bw = [0 0; 0 0; 1 0; 0 1]
R1 = Bw*I(2)*Bw'
R2 = 0.0001I(ny)
K = kalman(sys, R1, R2)
```

With our feedback gain `L` and the Kalman gain `K`, we form the controller using [`observer_controller`](@ref)
```@example PENDCART
controller = observer_controller(sys, L, K)
@assert isstable(controller)
@assert isstable(feedback(sys * controller))
```

We may have a look at the Nyquist plot and the gang-of-four to assess robustness margins. In this case we look at the loop transfer function at the input simply because this function is SISO while the standard output-loop transfer is MIMO. This will allow us to asses robustness w.r.t. input perturbations only
```@example PENDCART
nyquistplot(controller*sys, w, Ms_circles=[2.7], Mt_circles=[3], xlims=(-2, 2), ylims=(-1, 3))
```
The Nyquist plot shows a rather weak robustness margin, with a peak in the input sensitivity of about
```@example PENDCART
round(hinfnorm2(input_sensitivity(sys, controller))[1], digits=2) # hide
```
and a peak in the complementary sensitivity function of around
```@example PENDCART
round(hinfnorm2(input_comp_sensitivity(sys, controller))[1], digits=2) # hide
```
These can be verified by calling [`hinfnorm2`](@ref)
```@example PENDCART
hinfnorm2(input_comp_sensitivity(sys, controller))
```

!!! note "Hover information"
    If you plot with the Plotly backend, activated by calling `plotly()` if you have Plotly.jl installed, you can hover the mouse over the Nyquist curve and the gain circles to see frequency information etc. This is not possible when using the default GR backend, used in this documentation.

Also the gang-of-four indicate rather poor margins:
```@example PENDCART
gangoffourplot(sys, controller, w, xlabel="", sigma=false, titlefont=8)
```

### Robustification using Glover-McFarlane

In an attempt at improving this initial design, we call [`glover_mcfarlane`](@ref). This can be seen as a semi-automatic approach to robustifying an initial design, and will yield us an updated controller `Kgmf` with, hopefully, improved robustness properties.
```@example PENDCART
Kgmf, γ, info = glover_mcfarlane(sys, 1.05; W1=controller)
@assert isstable(Kgmf)
γ
```
The γ is an indication of the achieved robustness. A value of $γ < 4$ is typically desired, this time we did not quite achieve that, but will nevertheless proceed and look deeper into the robustness using other means.


!!! note "Controller order reduction"
    The Glover-McFarlane procedure often leads to high-order controllers. These controllers can sometimes be simplified by calling [`controller_reduction`](@ref), i.e., like this
    ```@example PENDCART
    Kgmfr, hs = controller_reduction(ExtendedStateSpace(sys), -Kgmf, 9) # Expects positive-feedback controller
    Kgmfr = -Kgmfr # Flip sign again
    Kgmfr.D .*= 0.0 # a hack to get better rolloff after reduction
    nothing # hide
    ```

### Robustness verification
We will now verify these designs in a number of ways. We start by inspecting sensitivity functions at the input, this function tells you how a load disturbance at the plant input translates to total plant input (including control action)
```@example PENDCART
f1 = bodeplot([controller*sys, Kgmf*sys], w, plot_title="Input Loop transfers", lab=["Pole placement" "" "GMF" ""]); vline!([2*4.85], sp=1, lab="Fundamental limitation", l=(:dash, :black))
f2 = nyquistplot([controller*sys, Kgmf*sys], xlims=(-4, 4), ylims=(-1, 5), Ms_circles=[2.7], Mt_circles=[3], lab=["Pole placement" "GMF"])
f3 = bodeplot(controller, w, lab="Pole placement")
bodeplot!(Kgmf, w, plot_title="Controllers", lab="GMF", legend=:bottomleft)
f4 = sigmaplot([
    input_sensitivity(sys, controller),
    input_sensitivity(sys, Kgmf)
    ], w, title="Input S", lab=["Pole placement" "GMF"], legend=:bottomright)
plot(f1,f2,f3,f4, size=(1000,1000))
```
We see that the Glover-Mcfarlane method increased the gain crossover frequency $\omega_{gc}$ slightly compared to the initial pole-placement controller, as well as lifted the phase even further. It also increased the roll-off, providing better filtering of high-frequency noise. However, it uses quite a bit more gain from the measurement of the pendulum angle.

The gang-of-four, shown below, looks slightly better for the GMF controller, shown in orange.
```@example PENDCART
gangoffourplot(sys, [controller, Kgmf], w, xlabel="", sigma=false, titlefontsize=8)
```
The robustified controller has better disturbance rejection ($P/(I + PC)$) and slightly lower peaks in the sensitivity and complementary sensitivity functions.

Inspecting the singular values of the output sensitivity, we see that the GMF controller reduces the peak and improves the disturbance rejection for the lower singular value, while leaving the upper singular value more or less where it is for low frequencies.
```@example PENDCART
sigmaplot(sensitivity.(Ref(sys), [controller, Kgmf]), w, lab=["Pole placement" "GMF"], legend=:bottomright)
```

Further verification of robustness properties can be conducted by inspecting the diskmargins at inputs and at outputs
```@example PENDCART
dmf1 = plot(diskmargin(sys*controller), title="Simultaneous Output diskmargin", lab="Pole placement")
dmf2 = plot(diskmargin(controller*sys), title="Input diskmargin", lab="Pole placement")
plot!(dmf1, diskmargin(sys*Kgmf), title="Simultaneous Output diskmargin", lab="GMF")
plot!(dmf2, diskmargin(Kgmf*sys), title="Input diskmargin", lab="GMF")
plot(dmf1, dmf2)
```
With the robustified controller, ee can tolerate a gain variation of about 1.6 at the plant input, but only 1.23 at the plant output. Please note that simultaneous margins can be quite conservative, it's much less likely that both outputs have equally large gain errors at the same time. One can also investigate the margins for one loop at a time using [`loop_diskmargin`](@ref).

### Simulation

Finally, it's time to simulate the system. First we simulate the response to a reference step for the cart position:
```@example PENDCART
plot([
    step(feedback(sys*controller)[:, 1], 8),
    step(feedback(sys*Kgmf)[:, 1], 8),
], ylab=["Pos" "Angle"], plot_title="Position command step response", lab=["Pole placement" "" "GMF" ""], legend=:bottomright)
```

Then we simulate the response to an impulsive disturbance acting on the cart (i.e., someone hit it with a hammer)
```@example PENDCART
plot([
    impulse(feedback(sys, controller), 8),
    impulse(feedback(sys, Kgmf), 8),
], ylab=["Pos" "Angle"], plot_title="Disturbance step response", lab=["Pole placement" "" "GMF" ""], legend=:bottomright)
```
This time, the controllers control also the cart position while keeping the pendulum stabilized. 


## Conclusion
We started out designing a PID controller and used the Bode plot to guide the tuning. While we ended up with a controller with good robustness margins, we had completely forgotten about the cart position and the controller turned out to not stabilize this "hidden state". We include this example here as an example of following a mostly sound procedure, leading to a robust controller, but failing to meet real-world constraints due to lack of observability.

The loop-shaping procedure yielded a controller that stabilized all states of the plant, but with questionable robustness margins. In practice, pole placement can be rather difficult and it's not always obvious where to place the poles to achieve a robust design. In this case, a robust design is very hard to achieve with a pole-placement controller without model augmentation, the poor robustness of the pole-placement controller compared to the PID controller is due to the low gain at low frequencies, indeed, the pole placement controller lacks integral action! See [Disturbance modeling and rejection with LQG controllers](@ref) for a tutorial on how to add integral action to state-feedback controllers by augmenting the system model with a disturbance model.

We looked at several different ways of quantifying robustness of a system with multiple outputs, and tried our luck with a procedure for automatic robustification, [`glover_mcfarlane`](@ref). In this case, the procedure worked and we got a slightly more robust controller as a result, this controller also increased the gain for low frequencies significantly, further indicating that the low-frequency gain was a source of problems for the pole-placement controller. The result of the Glover-McFarlane procedure may either be used directly as the final controller, or to provide insight into how the procedure modifies the existing controller in order to improve robustness.