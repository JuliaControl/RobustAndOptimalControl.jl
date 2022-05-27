# RobustAndOptimalControl.jl

This package is an extension to [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) that provides methods for robust and optimal analysis and synthesis of linear control systems. Some highlights:

- Named statespace systems ([`named_ss`](@ref)) where states, inputs and outputs are accessible by names rather than indices. This also facilitates creating complicated feedback interconnections using [`connect`](@ref).
- An interface to [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl). Call [`dss`](@ref) on a statespace system to get a descriptor system. We also forward some methods to implementations in DescriptorSystems.
- Robust/optimal design methods such as $H_{\infty}$, $H_{2}$, LQG and Glover-McFarlane.
- Robustness-related metrics such as [`nugap`](@ref) ($\nu$-gap), [`ncfmargin`](@ref), [`diskmargin`](@ref) etc.
- Uncertainty modeling with the $M\Delta$ framework (and more). Analsysis methods for this framework are still limited.
- Model augmentation.
- An [`ExtendedStateSpace`](@ref) type that represents a partitioned statespace system $w,u \rightarrow z,y$.



## Installation
```julia
pkg> add RobustAndOptimalControl
```

## Named systems
See [complicated-feedback example](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/complicated_feedback.jl)
- [`named_ss`](@ref)

Named systems can be indexed with their names, e.g.,
```julia
G[:y2, :u4]
```
but also using incomplete names, e.g., if `G` contains outputs `:y1, :y2, :y3, :z1, :z2`, the following retrieves the three outputs that has the prefix `:y`
```julia
G[:y, :] # Prefix matching is used if no exact match is found.
```

## Connecting systems together
See [complicated-feedback example](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/complicated_feedback.jl)
- [`connect`](@ref)

### Example
The following complicated feedback interconnection

```
                 yF
              ┌────────────────────────────────┐
              │                                │
    ┌───────┐ │  ┌───────┐ yR   ┌─────────┐    │    ┌───────┐
uF  │       │ │  │       ├──────►         │ yC │  uP│       │    yP
────►   F   ├─┴──►   R   │      │    C    ├────+────►   P   ├────┬────►
    │       │    │       │   ┌──►         │         │       │    │
    └───────┘    └───────┘   │  └─────────┘         └───────┘    │
                             │                                   │
                             └───────────────────────────────────┘
```
can be created by
```julia
F = named_ss(ssrand(1, 1, 2, proper=true), x=:xF, u=:uF, y=:yF)
R = named_ss(ssrand(1, 1, 2, proper=true), x=:xR, u=:uR, y=:yR)
C = named_ss(ssrand(1, 1, 2, proper=true), x=:xC, u=:uC, y=:yC)
P = named_ss(ssrand(1, 1, 3, proper=true), x=:xP, u=:uP, y=:yP)

addP = sumblock("uP = yF + yC") # Sum node before P
addC = sumblock("uC = yR - yP") # Sum node before C (drawn as two arrows into C in the diagram)

connections = [
    :yP => :yP # Output to input
    :uP => :uP # addP's output is called the same as P's input
    :yC => :yC
    :yF => :yF
    :yF => :uR
    :uC => :uC
    :yR => :yR
]
w1 = [:uF] # External inputs

G = connect([F, R, C, P, addP, addC], connections; w1)
```

If an external input is to be connected to multiple points, use a `splitter` to split up the signal into a set of unique names which are then used in the connections.



## Model augmentation
Add disturbance and performance models to your system model.

- [`add_disturbance`](@ref)
- [`add_measurement_disturbance`](@ref)
- [`add_input_differentiator`](@ref)
- [`add_output_differentiator`](@ref)
- [`add_input_integrator`](@ref)
- [`add_output_integrator`](@ref)
- [`add_low_frequency_disturbance`](@ref)
- [`add_resonant_disturbance`](@ref)

## $H_\infty$ and $H_2$ design
Examples are available in the [example folder](https://github.com/JuliaControl/RobustAndOptimalControl.jl/tree/master/examples).

- [`hinfsynthesize`](@ref)
- [`h2synthesize`](@ref)
- [`specificationplot`](@ref)
- [`glover_mcfarlane`](@ref)
- [`glover_mcfarlane_2dof`](@ref)
- [`hanus`](@ref)

### Example: Glover McFarlane design
This example will design a robust controller using the Glover-McFarlane method. This method requires the user to perform an initial loop-shaping design, i.e., by tuning a standard PI controller etc. The [`glover_mcfarlane`](@ref) method then takes the loop-shaping controller and the plant model and returns a robustified controller. This is example 9.3 from Skogestad, "Multivariable Feedback Control: Analysis and Design".
```@example GMF
using RobustAndOptimalControl, ControlSystems, Plots, Test
G = tf(200, [10, 1])*tf(1, [0.05, 1])^2     |> ss # Plant model
Gd = tf(100, [10, 1])                       |> ss # Disturbance model
W1 = tf([1, 2], [1, 1e-6])                  |> ss # Loop-shaping controller
K, γ, info = glover_mcfarlane(G, 1.1; W1)         # K is robustified controller
@test info.γmin ≈ 2.34 atol=0.005
Gcl = extended_gangoffour(G, K) # Form closed-loop system

fig1 = bodeplot([G, info.Gs, G*K], lab=["G" "" "Initial GK" "" "Robustified GK"])
fig2 = bodeplot(Gcl, lab=["S" "KS" "PS" "T"], plotphase=false) # Plot gang of four

# Simulate the response to a disturbance (Gd*feedback(1, G*K) = Gd*S is the closed-loop transfer function from an additive output disturbance)
fig3 = plot(step(Gd*feedback(1, info.Gs), 3), lab="Initial controller")
plot!(step(Gd*feedback(1, G*K), 3), lab="Robustified controller")
fig4 = nyquistplot([info.Gs, G*K], ylims=(-2,1), xlims=(-2, 1),
    Ms_circles = 1.5,
    lab = ["Initial controller" "Robustified controller"],
    title = "Loop transfers with and without robustified controller"
)
plot(fig1, fig2, fig3, fig4, titlefontsize=9, labelfontsize=9, size=(800, 640))
```

#### Example of controller reduction:
The order of the controller designed above can be reduced maintaining at least 2/3 of the robustness margin like this
```@example GMF
e,_ = ncfmargin(info.Gs, info.Ks)
Kr, hs, infor = baltrunc_coprime(info.Ks, n=info.Ks.nx)
n = findlast(RobustAndOptimalControl.error_bound(hs) .> 2e/3) # 2/3 e sets the robustness margin
Ksr, hs, infor = baltrunc_coprime(info.Ks; n)
@test ncfmargin(info.Gs, Ksr)[1] >= 2/3 * e
Kr = W1*Ksr
bodeplot([G*K, G*Kr], lab=["L original" "" "L Reduced" ""])
```
This gives a final controller `Kr` of order 2 instead of order 5, but a very similar robustness margin. You may also call
```@example GMF
controller_reduction_plot(info.Gs, info.Ks)
```
to help you select the controller order.


### Example: Glover McFarlane 2-dof design
In this example, we design a 2 degree-of-freedom controller using the Glover McFarlane method. This design method requires you to specify both a loop-shaping controller as well as a reference model. It's usually a good idea to let the reference model have the same number of poles as the system that is being controlled in order not not differentiate the references and introduce non-robustness.
```@example
using RobustAndOptimalControl, ControlSystems, Plots
P = tf([1, 5], [1, 2, 10]) # Plant
W1 = tf(1,[1, 0]) |> ss    # Loop shaping controller

Tref = tf(1, [1, 1])^2 |> ss # Reference model of same order as P

K1dof, γ1, info1 = glover_mcfarlane(ss(P), 1.1; W1)
K2dof, γ2, info2 = glover_mcfarlane_2dof(ss(P), Tref, 1.1, 1.1; W1)

G1 = feedback(P*K1dof)
G2 = info2.Gcl

w = exp10.(LinRange(-2, 2, 200))
fig1 = bodeplot(info2.K1, w, lab="Feedforward filter")
fig2 = plot([step(G1, 15), step(G2, 15), step(Tref, 15)], lab=["1-DOF" "2-DOF" "Tref"])
plot(fig1, fig2)
```

## LQG design
The main functionality for LQG design is exposed through [`LQGProblem`](@ref). See the docstring for an example.

## System analysis
- [`hinfnorm2`](@ref)
- [`linfnorm2`](@ref)
- [`hankelnorm`](@ref)
- [`h2norm`](@ref)
- [`nugap`](@ref)
- [`ncfmargin`](@ref)
- [`robstab`](@ref)
- [`ispassive`](@ref)
- [`passivity_index`](@ref)
- [`passivityplot`](@ref)



See also [Structured singular value and diskmargin](@ref) below
## Structured singular value and diskmargin
- [`structured_singular_value`](@ref). Note, this only handles diagonal complex perturbations at the moment.
- [`muplot`](@ref)
- [`diskmargin`](@ref)
- [`loop_diskmargin`](@ref)
- [`sim_diskmargin`](@ref)
- [`loop_scale`](@ref)
- [`loop_scaling`](@ref)

### Diskmargin example
The diskmargin can be visualized in several ways, as a region of allowed simultaneous gain and pahse variations:
```@example diskmargin
using RobustAndOptimalControl, ControlSystems, Plots
L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
plot(dm) # Plot the disk margin to illustrate maximum allowed simultaneous gain and phase variations.
```

As a Nyquist exclusion disk:
```@example diskmargin
nyquistplot(L)
plot!(dm, nyquist=true) # plot a nyquist exclusion disk. The Nyquist curve will be tangent to this disk at `dm.ω0`
nyquistplot!(dm.f0*L, lab="perturbed") # If we perturb the system with the worst-case perturbation `f0`, the curve will pass through the critical point -1.
```

And as a frequency-dependent margin
```@example diskmargin
w = exp10.(LinRange(-2, 2, 500))
dms = diskmargin(L, 0, w)
plot(dms)
```


## Closed-loop analysis
- [`output_sensitivity`](@ref)
- [`output_comp_sensitivity`](@ref)
- [`input_sensitivity`](@ref)
- [`input_comp_sensitivity`](@ref)
- [`G_CS`](@ref)
- [`G_PS`](@ref)
- [`gangoffour`](@ref)
- [`extended_gangoffour`](@ref)
- [`ncfmargin`](@ref)

## Model reduction

- [`baltrunc2`](@ref)
- [`frequency_weighted_reduction`](@ref)
- [`stab_unstab`](@ref)
- [`baltrunc_unstab`](@ref)
- [`baltrunc_coprime`](@ref)
- [`controller_reduction`](@ref)
- [`error_bound`](@ref)
- [`controller_reduction_plot`](@ref)