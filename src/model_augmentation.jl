"""
    add_disturbance(sys::StateSpace, Ad::Matrix, Cd::Matrix)

See CCS pp. 144

# Arguments:
- `sys`: System to augment
- `Ad`: The dynamics of the disturbance
- `Cd`: How the disturbance states affect the states of `sys`. This matrix has the shape (sys.nx, size(Ad, 1))

See also [`add_low_frequency_disturbance`](@ref), [`add_resonant_disturbance`](@ref)
"""
function add_disturbance(sys::AbstractStateSpace, Ad::AbstractMatrix, Cd::AbstractMatrix)
    A,B,C,D = ControlSystemsBase.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Ae = [A Cd; zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; zeros(T, size(Ad, 1), nu)]
    Ce = [C zeros(T, ny, size(Ad, 1))]
    De = D
    ss(Ae,Be,Ce,De,sys.timeevol)
end

"""
    add_measurement_disturbance(sys::StateSpace{Continuous}, Ad::Matrix, Cd::Matrix)

Create the system
```
Ae = [A 0; 0 Ad]
Be = [B; 0]
Ce = [C Cd]
```
"""
function add_measurement_disturbance(sys::AbstractStateSpace, Ad::AbstractMatrix, Cd::AbstractMatrix)
    A,B,C,D = ControlSystemsBase.ssdata(sys)
    T = eltype(A)
    @unpack nx,nu,ny = sys
    Ae = [A zeros(T, nx, size(Ad, 1)); zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; zeros(T, size(Ad, 1), nu)]
    Ce = [C Cd]
    De = D
    ss(Ae,Be,Ce,De,sys.timeevol)
end

"""
    add_low_frequency_disturbance(sys::StateSpace, Ai::Integer; ϵ = 0)

A disturbance affecting only state `Ai`.
"""
function add_low_frequency_disturbance(sys::AbstractStateSpace, Ai::Integer; ϵ=0)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    1 ≤ Ai ≤ nx || throw(ArgumentError("Ai must be a valid state index"))
    Cd = zeros(nx, 1)
    Cd[Ai] = 1
    Ad = -ϵ*I(nu)
    isdiscrete(sys) && (Ad += I)
    add_disturbance(sys, Ad, Cd)
end

"""
    add_low_frequency_disturbance(sys::StateSpace; ϵ = 0, measurement = false)
    add_low_frequency_disturbance(sys::StateSpace, Cd; ϵ = 0, measurement = false)

Augment `sys` with a low-frequency (integrating if `ϵ=0`) disturbance model.
If an integrating input disturbance is used together with an observer, the controller will have integral action.

- `Cd`: If adding an input disturbance. this matrix indicates how the disturbance states affect the states of `sys`, and defaults to `sys.B`. If `measurement=true`, this matrix indicates how the disturbance states affect the outputs of `sys`, and defaults to `I(sys.ny)`.

# Arguments:
- `ϵ`: Move the integrator pole `ϵ` into the stable region.
- `measurement`: If true, the disturbance is a measurement disturbance, otherwise it's an input diturbance. 
"""
function add_low_frequency_disturbance(sys::AbstractStateSpace, Cd::Union{Nothing, AbstractMatrix} = nothing; ϵ=0, measurement=false)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    if measurement
        Ad = -ϵ*I(ny)
        isdiscrete(sys) && (Ad += I)
        Cd === nothing && (Cd = I(ny))
        add_measurement_disturbance(sys, Ad, Cd)
    else
        Cd === nothing && (Cd = sys.B)
        Ad = -ϵ*I(size(Cd, 2)) # We use the size of Cd here in case not all inputs are augmented
        isdiscrete(sys) && (Ad += I)
        add_disturbance(sys, Ad, Cd)
    end
end

"""
    add_resonant_disturbance(sys::StateSpace{Continuous}, ω, ζ, Ai::Int; measurement = false)

Augment `sys` with a resonant disturbance model.

# Arguments:
- `ω`: Frequency
- `ζ`: Relative damping.
- `Ai`: The affected state
- `measurement`: If true, the disturbace is acting on the output, this will cause the controller to have zeros at ω (roots of poly s² + 2ζωs + ω²). If false, the disturbance is acting on the input, this will cause the controller to have poles at ω (roots of poly s² + 2ζωs + ω²).
"""
function add_resonant_disturbance(sys::AbstractStateSpace, ω, ζ, Ai::Integer; measurement=false)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    if measurement
        1 ≤ Ai ≤ sys.ny || throw(ArgumentError("Ai must be a valid output index"))
        Cd = zeros(ny, 2)
        Cd[Ai, 1] = 1
    else
        1 ≤ Ai ≤ sys.nx || throw(ArgumentError("Ai must be a valid state index"))
        Cd = zeros(nx, 2)
        Cd[Ai, 1] = 1
    end
    Ad = [-ζ -ω; ω -ζ]
    if isdiscrete(sys)
        Ad = exp(Ad * sys.Ts)
    end
    measurement ? add_measurement_disturbance(sys, Ad, Cd) : add_disturbance(sys, Ad, Cd)
end

"""
    add_resonant_disturbance(sys::AbstractStateSpace, ω, ζ, Bd::AbstractArray)

- `Bd`: The disturbance input matrix.
"""
function add_resonant_disturbance(sys::AbstractStateSpace, ω, ζ, Bd::AbstractArray)
    Ad = [-ζ -float(ω); ω -ζ]
    if isdiscrete(sys)
        Ad .*= sys.Ts
        Ad = exp(Ad)
    end
    add_disturbance(sys, Ad, [Bd zeros(sys.nx)])
end

"""
    add_differentiator(sys::StateSpace{<:Discrete})

Augment the output of `sys` with the numerical difference (discrete-time derivative) of output, i.e.,
`y_aug = [y; (y-y_prev)/sys.Ts]`
To add both an integrator and a differentiator to a SISO system, use
```julia
Gd = add_output_integrator(add_output_differentiator(G), 1)
```
"""
function add_output_differentiator(sys::AbstractStateSpace{<: Discrete}, diffsys=sys)
    A,B,C,D = ssdata(diffsys)
    all(iszero, D) || throw(ArgumentError("Can't add a differentiator to a system with non-zero D matrix. The system would not be proper."))
    C = C ./ sys.Ts
    Cd = C*(A-I)
    Dd = C*B
    A,B,C,D = ssdata(sys)
    ss(A, B, [C; Cd], [D; Dd], sys.timeevol)
end

function ControlSystemsBase.tf(M::AbstractArray{TransferFunction{TE,ControlSystemsBase.SisoRational{T}}}) where {TE, T<:Number}
    all(ControlSystemsBase.issiso, M) || throw(ArgumentError("To make a MIMO system out of several MIMO systems is not yet supported"))
    matrix = first.(getproperty.(M, :matrix))
    TransferFunction{TE,ControlSystemsBase.SisoRational{T}}(matrix, M[1].timeevol)
end

"""
    add_output_integrator(sys::StateSpace{<:Discrete}, ind = 1; ϵ = 0)

Augment the output of `sys` with the integral of output at index `ind`, i.e., 
`y_aug = [y; ∫y[ind]]`
To add both an integrator and a differentiator to a SISO system, use
```julia
Gd = add_output_integrator(add_output_differentiator(G), 1)
```

Note: numerical integration is subject to numerical drift. If the output of the system corresponds to, e.g., a velocity reference and the integral to position reference, consider methods for mitigating this drift.
"""
function add_output_integrator(sys::AbstractStateSpace{<: Discrete}, ind=1; ϵ=0, neg=false)
    int = tf(1.0*sys.Ts, [1, -(1-ϵ)], sys.Ts)
    neg && (int = int*(-1))
    𝟏 = tf(1.0,sys.Ts)
    𝟎 = tf(0.0,sys.Ts)
    M = [i==j ? 𝟏 : 𝟎 for i = 1:sys.ny, j = 1:sys.ny]
    M = [M; permutedims([i ∈ ind ? int : 𝟎 for i = 1:sys.ny])]
    nx = sys.nx
    nr = length(ind)
    p = [(1:nx).+nr; 1:nr]
    T = (1:nx+nr) .== p'
    similarity_transform(tf(M)*sys, T)
end

function add_output_integrator(sys::AbstractStateSpace{Continuous}, ind=1; ϵ=0, neg=false)
    int = tf(1.0, [1, ϵ])
    neg && (int = int*(-1))
    𝟏 = tf(1.0)
    𝟎 = tf(0.0)
    M = [i==j ? 𝟏 : 𝟎 for i = 1:sys.ny, j = 1:sys.ny]
    M = [M; permutedims([i ∈ ind ? int : 𝟎 for i = 1:sys.ny])]
    nx = sys.nx
    nr = length(ind)
    p = [(1:nx).+nr; 1:nr]
    T = (1:nx+nr) .== p'
    similarity_transform(tf(M)*sys, T)
end

"""
    add_input_integrator(sys::StateSpace, ui = 1, ϵ = 0)

Augment the output of `sys` with the integral of input at index `ui`, i.e., 
`y_aug = [y; ∫u[ui]]`
See also [`add_low_frequency_disturbance`](@ref)
"""
function add_input_integrator(sys::AbstractStateSpace, ui=1; ϵ=0)
    A,B,C,D = ControlSystemsBase.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    1 ≤ ui ≤ nu || throw(ArgumentError("ui must be a valid input index"))
    Cd = zeros(T, 1, nx+1)
    Cd[end] = 1
    Bd = zeros(T, 1, nu)
    Bd[ui] = ControlSystemsBase.isdiscrete(sys) ? sys.Ts : 1
    Ad = -ϵ*I(1)
    isdiscrete(sys) && (Ad += I)

    Ae = [A zeros(T, nx, 1); zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; Bd]
    Ce = [[C zeros(T, ny, 1)]; Cd]
    De = [D; zeros(T, 1, nu)]
    ss(Ae,Be,Ce,De,sys.timeevol)

end


"""
    add_input_differentiator(sys::StateSpace, ui = 1:sys.nu; goodwin=false)

Augment the output of `sys` with the difference `u(k+1)-u(k)`

# Arguments:
- `ui`: An index or vector of indices indicating which inputs to differentiate.
- `goodwin`: If true, the difference operator will use the Goodwin δ operator, i.e., `(u(k+1)-u(k)) / sys.Ts`.

The augmented system will have the matrices
```
[A 0; 0 0]  [B; I]  [C 0; 0 -I]  [D; I]
```
with `length(ui)` added states and outputs.
"""
function add_input_differentiator(sys::AbstractStateSpace{<:Discrete}, ui=1:sys.nu; goodwin=false)
    A,B,C,D = ControlSystemsBase.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    all(1 .≤ ui .≤ nu) || throw(ArgumentError("ui must be a valid input index"))
    nnu = length(ui) # number of new states and outputs

    den = goodwin ? 1/sys.Ts : 1

    Cd = zeros(T, nnu, nx+nnu)
    Cd[:, nx+1:end] .= -den*I(nnu)
    Bd = zeros(T, nnu, nu)
    for (i, ui) in enumerate(ui)
        Bd[i, ui] = 1
    end
    Ad = zeros(nnu, nnu)
    Dd = den*I(nnu)

    Ae = [A zeros(T, nx, nnu); zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; Bd]
    Ce = [[C zeros(T, ny, nnu)]; Cd]
    De = [D; Dd]
    ss(Ae,Be,Ce,De,sys.timeevol)

end

# using ControlSystemsBase.DemoSystems
# sys = DemoSystems.resonant()
# sys2 = add_low_frequency_disturbance(sys, 2)
# sys25 = add_low_frequency_disturbance(sys)
# sys3 = add_resonant_disturbance(sys, 1, 0.5, 1)
# ss([0], [0], [1], 1)*sys
# sys + ss([0.0], [0], [1], 1)