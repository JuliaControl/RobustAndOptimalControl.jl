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

# Named systems
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

# Connecting systems together
See [complicated-feedback example](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/complicated_feedback.jl)
- [`connect`](@ref)

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



# Model augmentation
TODO.

Add disturbance and performance models to your system model.

- [`add_disturbance`](@ref)
- [`add_measurement_disturbance`](@ref)
- [`add_input_differentiator`](@ref)
- [`add_output_differentiator`](@ref)
- [`add_input_integrator`](@ref)
- [`add_output_integrator`](@ref)
- [`add_low_frequency_disturbance`](@ref)
- [`add_resonant_disturbance`](@ref)

# $H_\infty$ and $H_2$ design
TODO
Examples are available in the [example folder](https://github.com/JuliaControl/RobustAndOptimalControl.jl/tree/master/examples).

- [`hinfsynthesize`](@ref)
- [`h2synthesize`](@ref)
- [`glover_mcfarlane`](@ref)
- [`glover_mcfarlane_2dof`](@ref)
# LQG design
TODO
- [`LQGProblem`](@ref)
# Structured singular value and diskmargin
- [`structured_singular_value`](@ref). Note, this only handles diagonal complex perturbations at the moment.
- [`diskmargin`](@ref)
- [`loop_diskmargin`](@ref)
- [`sim_diskmargin`](@ref)


# Closed-loop analysis
- [`output_sensitivity`](@ref)
- [`output_comp_sensitivity`](@ref)
- [`input_sensitivity`](@ref)
- [`input_comp_sensitivity`](@ref)
- [`G_CS`](@ref)
- [`G_PS`](@ref)
- [`gangoffour`](@ref)

