# RobustAndOptimalControl.jl

## Installation
```julia
pkg> add RobustAndOptimalControl
```

# Named systems
[`named_ss`](@ref)

# Connecting systems together
[`connect`](@ref)

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
addC = sumblock("uC = yR - yP") # Sum node before C

connections = [
    :yP => :yP # Output to input
    :uP => :uP
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
# Uncertainty modeling
TODO: add examples

- [`neglected_delay`](@ref)
- [`neglected_lag`](@ref)
- [`gain_and_delay_uncertainty`](@ref)
- [`makeweight`](@ref)
- [`fit_complex_perturbations`](@ref)

# Model augmentation
TODO

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
# LQG design
TODO
- [`LQGProblem`](@ref)
# Structured singular value and diskmargin
TODO

# Closed-loop analysis
- [`output_sensitivity`](@ref)
- [`output_comp_sensitivity`](@ref)
- [`input_sensitivity`](@ref)
- [`input_comp_sensitivity`](@ref)
- [`G_CS`](@ref)
- [`G_PS`](@ref)
- [`gangoffour`](@ref)

# Exported functions and types
## Index

```@index
```
```@autodocs
Modules = [RobustAndOptimalControl]
Private = false
```
