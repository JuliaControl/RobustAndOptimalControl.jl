"""
    add_disturbance(sys::AbstractStateSpace, Ad::AbstractMatrix, Cd::AbstractMatrix)

See CCS pp. 144

# Arguments:
- `sys`: System to augment
- `Ad`: The dynamics of the disturbance
- `Cd`: How the disturbance states affect the states of `sys`. This matrix has the shape (sys.nx, size(Ad, 1))

See also `add_low_frequency_disturbance, add_resonant_disturbance`
"""
function add_disturbance(sys::AbstractStateSpace, Ad::AbstractMatrix, Cd::AbstractMatrix)
    A,B,C,D = ControlSystems.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Ae = [A Cd; zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; zeros(T, size(Ad, 1), nu)]
    Ce = [C zeros(T, ny, size(Ad, 1))]
    De = D
    ss(Ae,Be,Ce,De,sys.timeevol)
end

function add_measurement_disturbance(sys::AbstractStateSpace{Continuous}, Ad::AbstractMatrix, Cd::AbstractMatrix)
    A,B,C,D = ControlSystems.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Ae = [A zeros(T, nx, size(Ad, 1)); zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; zeros(T, size(Ad, 1), nu)]
    Ce = [C Cd]
    De = D
    ss(Ae,Be,Ce,De)
end

function add_low_frequency_disturbance(sys::AbstractStateSpace, Ai::Integer; ϵ=0)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    1 ≤ Ai ≤ nx || throw(ArgumentError("Ai must be a valid state index"))
    Cd = zeros(nx, 1)
    Cd[Ai] = 1
    Ad = -ϵ*I(nu)
    isdiscrete(sys) && (Ad += I)
    add_disturbance(sys, Ad, Cd)
end

function add_low_frequency_disturbance(sys::AbstractStateSpace; ϵ=0, measurement=false)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    Ad = -ϵ*I(nu)
    isdiscrete(sys) && (Ad += I)
    if measurement
        Cd = I(nu)
        add_measurement_disturbance(sys, Ad, Cd)
    else
        Cd = sys.B
        add_disturbance(sys, Ad, Cd)
    end
end

function add_resonant_disturbance(sys::AbstractStateSpace{Continuous}, ω, ζ, Ai::Integer; measurement=false)
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
    measurement ? add_measurement_disturbance(sys, Ad, Cd) : add_disturbance(sys, Ad, Cd)
end

"""
    add_differentiator(sys::AbstractStateSpace{<:Discrete})

Augment the output of `sys` with the numerical difference (discrete-time derivative) of output, i.e.,
`y_aug = [y; (y-y_prev)/sys.Ts]`
To add both an integrator and a differentiator to a SISO system, use
```

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

function ControlSystems.tf(M::AbstractArray{TransferFunction{TE,ControlSystems.SisoRational{T}}}) where {TE, T<:Number}
    all(ControlSystems.issiso, M) || throw(ArgumentError("To make a MIMO system out of several MIMO systems is not yet supported"))
    matrix = first.(getproperty.(M, :matrix))
    TransferFunction{TE,ControlSystems.SisoRational{T}}(matrix, M[1].timeevol)
end

"""
add_output_integrator(sys::AbstractStateSpace{<:Discrete}, ind = 1; ϵ = 0)

Augment the output of `sys` with the integral of output at index `ind`, i.e., 
`y_aug = [y; ∫y[ind]]`
To add both an integrator and a differentiator to a SISO system, use
```
Gd = add_output_integrator(add_output_differentiator(G), 1)
```

Note: numerical integration is subject to numerical drift. If the output of the system corresponds to, e.g., a velocity reference and the integral to position reference, consider methods for mitigating this drift.
"""
function add_output_integrator(sys::AbstractStateSpace{<: Discrete}, ind=1; ϵ=0)
    int = tf(1.0, [1, -(1-ϵ)], sys.Ts)
    𝟏 = tf(1.0,sys.Ts)
    𝟎 = tf(0.0,sys.Ts)
    M = [i==j ? 𝟏 : 𝟎 for i = 1:sys.ny, j = 1:sys.ny]
    M = [M; permutedims([i==ind ? int : 𝟎 for i = 1:sys.ny])]
    tf(M)*sys
end

"""
    add_input_integrator(sys::AbstractStateSpace, ui = 1, ϵ = 0)

Augment the output of `sys` with the integral of input at index `ui`, i.e., 
`y_aug = [y; ∫u[ui]]`
"""
function add_input_integrator(sys::AbstractStateSpace, ui=1; ϵ=0)
    A,B,C,D = ControlSystems.ssdata(sys)
    T = eltype(A)
    nx,nu,ny = sys.nx,sys.nu,sys.ny
    1 ≤ ui ≤ nu || throw(ArgumentError("ui must be a valid input index"))
    Cd = zeros(T, 1, nx+1)
    Cd[end] = 1
    Bd = zeros(T, 1, nu)
    Bd[ui] = 1
    Ad = -ϵ*I(1)
    isdiscrete(sys) && (Ad += I)

    Ae = [A zeros(T, nx, 1); zeros(T, size(Ad, 1), nx) Ad]
    Be = [B; Bd]
    Ce = [[C zeros(T, ny, 1)]; Cd]
    De = [D; zeros(T, 1, nu)]
    ss(Ae,Be,Ce,De,sys.timeevol)

end


# using ControlSystems.DemoSystems
# sys = DemoSystems.resonant()
# sys2 = add_low_frequency_disturbance(sys, 2)
# sys25 = add_low_frequency_disturbance(sys)
# sys3 = add_resonant_disturbance(sys, 1, 0.5, 1)
# ss([0], [0], [1], 1)*sys
# sys + ss([0.0], [0], [1], 1)