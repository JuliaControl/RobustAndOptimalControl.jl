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

Augment `sys` with a single low-frequency (integrating if `ϵ=0`) disturbance state
affecting state index `Ai` of `sys`.

# Arguments:
- `Ai`: Index of the plant state the disturbance is added to. Must satisfy `1 ≤ Ai ≤ sys.nx`.
- `ϵ`: Move the integrator pole `ϵ` into the stable region (continuous: pole at `-ϵ`; discrete: pole at `1-ϵ`).
"""
function add_low_frequency_disturbance(sys::AbstractStateSpace, Ai::Integer; ϵ=0)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    1 ≤ Ai ≤ nx || throw(ArgumentError("Ai must be a valid state index"))
    Cd = zeros(nx, 1)
    Cd[Ai] = 1
    Ad = fill(-float(ϵ), 1, 1)
    isdiscrete(sys) && (Ad .+= 1)
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

Augment `sys` with a resonant disturbance model. The added disturbance dynamics have
eigenvalues `-ζ ± iω`, i.e. `ζ` is the decay rate and `ω` is the damped (oscillation)
frequency. The characteristic polynomial of the continuous-time disturbance dynamics is
`s² + 2ζs + (ζ² + ω²)`. For an undamped oscillator at frequency `ω`, choose `ζ = 0`.

# Arguments:
- `ω`: Damped (oscillation) frequency, i.e. the imaginary part of the disturbance poles.
- `ζ`: Decay rate, i.e. the real part of the disturbance poles (units: inverse time).
- `Ai`: The affected state
- `measurement`: If true, the disturbance acts on the output, causing the controller to have zeros near the disturbance poles. If false, the disturbance acts on the input, causing the controller to have poles near the disturbance poles.
"""
function add_resonant_disturbance(sys::AbstractStateSpace, ω, ζ, Ai::Integer; measurement=false)
    A, _, _, _ = ControlSystemsBase.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    if measurement
        1 ≤ Ai ≤ sys.ny || throw(ArgumentError("Ai must be a valid output index"))
        Cd = zeros(T, ny, 2)
        Cd[Ai, 1] = 1
    else
        1 ≤ Ai ≤ sys.nx || throw(ArgumentError("Ai must be a valid state index"))
        Cd = zeros(T, nx, 2)
        Cd[Ai, 1] = 1
    end
    Ad = T[-ζ -ω; ω -ζ]
    if isdiscrete(sys)
        Ad = exp(Ad * sys.Ts)
    end
    measurement ? add_measurement_disturbance(sys, Ad, Cd) : add_disturbance(sys, Ad, Cd)
end

"""
    add_resonant_disturbance(sys::AbstractStateSpace, ω, ζ, Bd::AbstractArray; measurement = false)

Augment `sys` with a resonant disturbance whose injection into `sys` is described by `Bd`.
See the integer-`Ai` method for the meaning of `ω` and `ζ`.

# Arguments:
- `Bd`: Disturbance injection matrix.
    - If `measurement = false`, `Bd` indicates how the disturbance states affect the plant
      states. It must have `sys.nx` rows and either `1` or `2` columns. With one column,
      only the first (cosine-like) resonant state injects into the plant; with two columns,
      both resonant states inject.
    - If `measurement = true`, `Bd` indicates how the disturbance states affect the plant
      outputs and must have `sys.ny` rows and `2` columns.
"""
function add_resonant_disturbance(sys::AbstractStateSpace, ω, ζ, Bd::AbstractArray; measurement=false)
    A, _, _, _ = ControlSystemsBase.ssdata(sys)
    T = eltype(A)
    Ad = T[-ζ -ω; ω -ζ]
    if isdiscrete(sys)
        Ad = exp(Ad * sys.Ts)
    end
    if measurement
        size(Bd, 1) == sys.ny || throw(ArgumentError("Bd must have sys.ny=$(sys.ny) rows in the measurement case, got $(size(Bd, 1))"))
        size(Bd, 2) == 2 || throw(ArgumentError("Bd must have 2 columns in the measurement case (one per resonant state), got $(size(Bd, 2))"))
        add_measurement_disturbance(sys, Ad, Bd)
    else
        size(Bd, 1) == sys.nx || throw(ArgumentError("Bd must have sys.nx=$(sys.nx) rows, got $(size(Bd, 1))"))
        nc = size(Bd, 2)
        nc == 1 || nc == 2 || throw(ArgumentError("Bd must have 1 or 2 columns, got $nc"))
        Cd = nc == 2 ? Bd : [Bd zeros(T, sys.nx)]
        add_disturbance(sys, Ad, Cd)
    end
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
    add_output_integrator(sys::StateSpace, ind = 1; ϵ = 0, neg = false)

Augment the output of `sys` with the integral of the outputs at indices `ind`, i.e.,
`y_aug = [y; ∫y[ind]]`. One integrator state is added per entry of `ind`, in the order the
indices are given, and the integrator states are appended after the states of `sys`:
```math
\\begin{bmatrix} \\dot{x} \\\\ \\dot{x_i} \\end{bmatrix} =
\\begin{bmatrix} A & 0 \\\\ C_i & -ϵI \\end{bmatrix}
\\begin{bmatrix} x \\\\ x_i \\end{bmatrix} +
\\begin{bmatrix} B \\\\ D_i \\end{bmatrix} u
```
where `Cᵢ = C[ind, :]` and `Dᵢ = D[ind, :]`. For a discrete-time system, the integrator states
obey `xᵢ(k+1) = (1-ϵ)xᵢ(k) + Ts*y[ind](k)` (forward Euler), so that `xᵢ` approximates the
time integral of `y[ind]` in both time domains.

To add both an integrator and a differentiator to a SISO system, use
```julia
Gd = add_output_integrator(add_output_differentiator(G), 1)
```

# Arguments:
- `ind`: Output indices to integrate. Accepts an `Integer`, an `AbstractVector{<:Integer}` or an `AbstractRange`.
- `ϵ`: Move the integrator poles into the stable region, to `-ϵ` in continuous time and to `1-ϵ` in discrete time.
- `neg`: Negate the added outputs, i.e., `y_aug = [y; -∫y[ind]]`. This affects the added output rows only, never the integrator state dynamics.

Note: numerical integration is subject to numerical drift. If the output of the system corresponds to, e.g., a velocity reference and the integral to position reference, consider methods for mitigating this drift.
"""
function add_output_integrator(sys::AbstractStateSpace, ind=1; ϵ=0, neg=false)
    inds = ind isa Integer ? (ind:ind) : ind
    all(i -> 1 ≤ i ≤ sys.ny, inds) || throw(ArgumentError("All output indices in ind = $ind must be in 1:$(sys.ny)"))
    A, B, C, D = ssdata(sys)
    nx, nu, ny = sys.nx, sys.nu, sys.ny
    nr = length(inds)
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D), typeof(ϵ), Float64)
    # The integrator state is a genuine time integral in both time domains, so that the
    # weights applied to it by, e.g., `lqi` carry the same meaning for a continuous-time
    # system and its discretization.
    h = isdiscrete(sys) ? T(sys.Ts) : one(T)
    λ = isdiscrete(sys) ? 1 - ϵ : -ϵ
    Aa = T[A zeros(nx, nr); h*C[inds, :] λ*I(nr)]
    Ba = T[B; h*D[inds, :]]
    Ci = neg ? -I(nr) : I(nr)
    Ca = T[C zeros(ny, nr); zeros(nr, nx) Ci]
    Da = T[D; zeros(nr, nu)]
    ss(Aa, Ba, Ca, Da, sys.timeevol)
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
    C_int_row = zeros(T, 1, nx+1)
    C_int_row[end] = 1
    B_int_row = zeros(T, 1, nu)
    B_int_row[ui] = ControlSystemsBase.isdiscrete(sys) ? sys.Ts : 1
    A_int = -ϵ*I(1)
    isdiscrete(sys) && (A_int += I)

    Ae = [A zeros(T, nx, 1); zeros(T, size(A_int, 1), nx) A_int]
    Be = [B; B_int_row]
    Ce = [[C zeros(T, ny, 1)]; C_int_row]
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