
import Base.getindex
import ControlSystems.numeric_type

"""
    G = LQGProblem(sys::ExtendedStateSpace, Q1, Q2, R1, R2; qQ=0, qR=0, SQ=nothing, SR=nothing)

Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)`
where the controller is of LQG-type. The controller is specified by weight matrices `Q1,Q2`
that penalizes state deviations and control signal variance respectively, and covariance
matrices `R1,R2` which specify state drift and measurement covariance respectively.

`sys` is an extended statespace object where the upper channel corresponds to disturbances to performance variables (w→z), and the lower channel corresponds to inputs to outputs (u→y), such that `lft(sys, K)` forms the closed-loop transfer function from external inputs/disturbances to performance variables. 

`qQ` and `qR` can be set to incorporate loop-transfer recovery, i.e.,
```julia
L = lqr(A, B, Q1+qQ*C'C, Q2)
K = kalman(A, C, R1+qR*B*B', R2)
```
Increasing `qQ` will add more cost in output direction, e.g., encouraging the use of cheap control, while
increasing `qR` adds fictious dynamics noise, makes the observer faster in the direction we control.

# Example

```julia
s = tf("s")
P = [1/(s+1) 2/(s+2); 1/(s+3) 1/(s-1)]
sys = ss(P)
eye(n) = Matrix{Float64}(I,n,n) # For convinience

qQ = 1
qR = 1
Q1 = 10eye(4)
Q2 = 1eye(2)
R1 = 1eye(6)
R2 = 1eye(2)

G = LQGProblem(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR)

Gcl = G.cl
T = G.T
S = G.S
sigmaplot([S,T],exp10.(range(-3, stop=3, length=1000)))
stepplot(Gcl)
```

"""
struct LQGProblem
    sys::ExtendedStateSpace
    Q1::AbstractMatrix
    Q2::AbstractMatrix
    R1::AbstractMatrix
    R2::AbstractMatrix
    qQ::Real
    qR::Real
    SQ::AbstractMatrix
    SR::AbstractMatrix
end

ControlSystems.isdiscrete(l::LQGProblem) = ControlSystems.isdiscrete(l.sys)
ControlSystems.iscontinuous(l::LQGProblem) = ControlSystems.iscontinuous(l.sys)

function LQGProblem(
    sys::AbstractStateSpace,
    Q1::AbstractVecOrMat,
    Q2::AbstractVecOrMat,
    R1::AbstractVecOrMat,
    R2::AbstractVecOrMat;
    qQ = 0,
    qR = 0,
    SQ = nothing,
    SR = nothing,
    kwargs...,
)
    sys isa ExtendedStateSpace || (sys = ExtendedStateSpace(sys, B1=I, C1=I))
    @unpack B1, C1, C2 = sys
    Q1 = Q1 isa AbstractVector ? diagm(Q1) : Q1
    Q2 = Q2 isa AbstractVector ? diagm(Q2) : Q2
    R1 = R1 isa AbstractVector ? diagm(R1) : R1
    R2 = R2 isa AbstractVector ? diagm(R2) : R2
    size(Q1, 1) == size(C1,1) || throw(ArgumentError("The size of Q1 is determined by C1, not by the state."))
    size(R2, 1) == size(C2,1) || throw(ArgumentError("The size of R2 is determined by C2, not by the state."))
    size(R1, 1) == size(B1,2) || throw(ArgumentError("The size of R1 is determined by B1, not by the state."))
    SQ === nothing && (SQ = zeros(size(sys.B2)))
    SR === nothing && (SR = zeros(size(sys.C2')))

    LQGProblem(sys, Q1, Q2, R1, R2, qQ, qR, SQ, SR; kwargs...)
end 

"""
    LQGProblem(P::ExtendedStateSpace)

If only an `ExtendedStateSpace` system is provided, e.g. from `hinfpartition`, the system `P` is assumed to correspond to the H₂ optimal control problem with
```
C1'C1    = Q1
D12'D12  = Q2
SQ       = C1'D12 # Cross term

B1*B1'   = R1
D21*D21' = R2
SR       = B1*D21' # Cross term
```
and an `LQGProblem` with the above covariance matrices is returned. The system
description in the returned LQGProblem will have `B1 = C1 = I`.
See Ch. 14 in Robust and optimal control for reference. 

# Example:
All the following ways of obtaining the H2 optimal controller are (almost) equivalent
```julia
using Test
G = ss(tf(1, [10, 1]))
WS = tf(1, [1, 1e-6]) 
WU = makeweight(1e-2, 0.1, 100) 
Gd = hinfpartition(G, WS, WU, [])

K, Gcl = h2synthesize(Gd)              # First option, using H2 formulas
K2, Gcl2 = h2synthesize(Gd, 1000)      # Second option, using H∞ formulas with large γ

lqg = LQGProblem(Gd)                   # Third option, convert to an LQGProblem and obtain controller
K3 = -observer_controller(lqg)

@test h2norm(lft(Gd, K )) ≈ 3.0568 atol=1e-3
@test h2norm(lft(Gd, K2)) ≈ 3.0568 atol=1e-3
@test h2norm(lft(Gd, K3)) ≈ 3.0568 atol=1e-3
```
"""
function LQGProblem(P::ExtendedStateSpace)
    @unpack A, B1, B2, C1, C2, D11, D12, D21, D22 = P
    # all(iszero, D22) || throw(ArgumentError("Non-zero D22 not handled."))
    all(iszero, D11) || throw(ArgumentError("Non-zero D11 not handled."))
    # all(iszero, D12'C1) || throw(ArgumentError("D12'C1 should be 0")) # One could get around this if using the cross-term in ARE, see 
    # all(iszero, D21*B1') || throw(ArgumentError("D21*B1' should be 0"))
    Q1 = C1'C1
    Q2 = D12'D12
    R1 = B1*B1'
    R2 = D21*D21'
    SR = B1*D21'
    SQ = C1'D12
    B1 = I(P.nx)
    C1 = I(P.nx)

    P = ss(A, B1, B2, C1, C2; D22, Ts = P.timeevol)
    # P = ss(A, B1, B2, C1, C2; Ts = P.timeevol)
    LQGProblem(P, Q1, Q2, R1, R2; SQ, SR)
end

function ControlSystems.kalman(l::LQGProblem)
    @unpack A, C2, B1, R1, qR, B2, R2, SR = l
    K = kalman(l.timeevol, A, C2, Hermitian(B1*R1*B1' + qR * B2 * B2'), R2, SR)
end

function ControlSystems.lqr(l::LQGProblem)
    @unpack A, B2, C1, Q1, qQ, C2, Q2, SQ = l 
    L = lqr(l.timeevol, A, B2, Hermitian(C1'Q1*C1 + qQ * C2'C2), Q2, SQ)
end

"""
    lqr3(P::AbstractStateSpace, Q1::AbstractMatrix, Q2::AbstractMatrix, Q3::AbstractMatrix)

Calculate the feedback gain of the discrete LQR cost function augmented with control differences
```math
x^{T} Q_1 x + u^{T} Q_2 u + Δu^{T} Q_3 Δu, \\quad
Δu = u(k) - u(k-1)
```
"""
function lqr3(P::AbstractStateSpace{<:Discrete}, Q1::AbstractMatrix, Q2::AbstractMatrix, Q3::AbstractMatrix)
    Pd = add_input_differentiator(P)
    S = zeros(Pd.nx, P.nu)
    S[P.nx+1:end, :] = -Q3
	X, _, L = MatrixEquations.ared(Pd.A, Pd.B, Q2+Q3, cat(Q1, Q3, dims=(1,2)), S) # ME has cost matrices reversed
    L[:, 1:P.nx]
end

"""
    dare3(P::AbstractStateSpace, Q1::AbstractMatrix, Q2::AbstractMatrix, Q3::AbstractMatrix; full=false)

Solve the discrete-time algebraic Riccati equation for a discrete LQR cost augmented with control differences
```math
x^{T} Q_1 x + u^{T} Q_2 u + Δu^{T} Q_3 Δu, \\quad
Δu = u(k) - u(k-1)
```

If `full`, the returned matrix will include the state `u(k-1)`, otherwise the returned matrix will be of the same size as `Q1`.
"""
function dare3(P::AbstractStateSpace{<:Discrete}, Q1::AbstractMatrix, Q2::AbstractMatrix, Q3::AbstractMatrix; full=false)
    # The reference cited in MatrixEquations.ared, W.F. Arnold, III and A.J. Laub, Generalized Eigenproblem Algorithms and Software for Algebraic Riccati Equations
    # defines the cost function as x'Q1x + u'Q2u + 2x'Su.
    # The Δu term expands to u+'Q2u + u'Q2u - 2u Q3 u+, so the factor 2 is already accounted for
    Pd = add_input_differentiator(P)
    S = zeros(Pd.nx, P.nu)
    S[P.nx+1:end, :] = -Q3
	X, _, L = MatrixEquations.ared(Pd.A, Pd.B, Q2+Q3, cat(Q1, Q3, dims=(1,2)), S) # ME has cost matrices reversed
    full ? X : X[1:P.nx, 1:P.nx]
end

dare3(A::AbstractMatrix, B, Q1, Q2, Q3::AbstractMatrix) = dare3(ss(A, B, I(size(A,1)), 0, 1), Q1, Q2, Q3)

function static_gain_compensation(l::LQGProblem, L = lqr(l))
    @unpack A, C1, B1, B2, D11, D12  = l
    pinv(D12 - (C1 - D12*L) * inv(A - B2*L) * B2)
end

function static_gain_compensation(A, B, C, D, L)
    pinv(D - (C - D*L) * inv(A - B*L) * B) # if D is nonzero
    # pinv(C * ((B * L - A) \ B))
end


"""
    extended_controller(K::AbstractStateSpace)

Takes a controller and returns an `ExtendedStateSpace` version which has augmented input `[r; y]` and output `y` (`z` output is 0-dim).
"""
function extended_controller(K::AbstractStateSpace)
    nx,nu,ny = K.nx, K.nu, K.ny
    A,B,C,D = ssdata(K)
    @error("This has not been verified")

    B1 = zeros(nx, nx) # dynamics not affected by r
    B2 = B # input y
    D21 = C#   K*r
    C2 = -C # - K*y
    C1 = zeros(0, nx)
    ss(A, B1, B2, C1, C2; D21, Ts = K.timeevol)
end

"""
    extended_controller(l::LQGProblem, L = lqr(l), K = kalman(l))

Returns an expression for the controller that is obtained when state-feedback `u = -L(xᵣ-x̂)` is combined with a Kalman filter with gain `K` that produces state estimates x̂. The controller is an instance of `ExtendedStateSpace` where `C2 = -L, D21 = L` and `B2 = K`.

The returned system has *inputs* `[xᵣ; y]` and outputs the control signal `u`. If a reference model `R` is used to generate state references `xᵣ`, the controller from `e = ry - y -> u` is given by
```julia
Ce = extended_controller(l)
Ce = named_ss(Ce; x = :xC, y = :u, u = [R.y; :y^l.ny]) # Name the inputs of Ce the same as the outputs of `R`.
connect([R, Ce]; u1 = R.y, y1 = R.y, w1 = [:ry^l.ny, :y^l.ny], z1=[:u])
```

Since the negative part of the feedback is built into the returned system, we have
```julia
C = observer_controller(l)
Ce = extended_controller(l)
system_mapping(Ce) == -C
```
"""
function extended_controller(l::LQGProblem, L = lqr(l), K = kalman(l))
    A,B,C,D = ssdata(system_mapping(l))
    Ac = A - B*L - K*C + K*D*L # 8.26b
    nx = l.nx
    B1 = zeros(nx, nx) # dynamics not affected by r
    B2 = K # input y
    D21 = L #   L*xᵣ # should be D21?
    C2 = -L # - L*x̂
    C1 = zeros(0, nx)
    ss(Ac, B1, B2, C1, C2; D21, Ts = l.timeevol)
end

"""
    observer_controller(l::LQGProblem, L = lqr(l), K = kalman(l))

Returns an expression for the feedback controller `u = Cy` that is obtained when state-feedback `u = -Lx̂` is combined with a Kalman filter with gain `K` that produces state estimates x̂.

Note: the transfer function returned is only a representation of the controller in the simple setting described above, e.g., it is not valid if the actual input contains anything that is not produced by a pure feedback from observed states. To obtain a controller that takes references into account, see `extended_controller`.
"""
function ControlSystems.observer_controller(l::LQGProblem, L::AbstractMatrix = lqr(l), K::AbstractMatrix = kalman(l))
    A,B,C,D = ssdata(system_mapping(l, identity))
    Ac = A - B*L - K*C + K*D*L # 8.26b
    Bc = K
    Cc = L
    Dc = 0
    ss(Ac, Bc, Cc, Dc, l.timeevol)
end

function ff_controller(l::LQGProblem, L = lqr(l), K = kalman(l))
    Ae,Be,Ce,De = ssdata(system_mapping(l, identity))
    Ac = Ae - Be*L - K*Ce + K*De*L # 8.26c
    Bc = Be*static_gain_compensation(l, L)
    Cc = L
    Dc = 0
    return 1 - ss(Ac, Bc, Cc, Dc, l.timeevol)
end

"""
    closedloop(l::LQGProblem, L = lqr(l), K = kalman(l))

Closed-loop system as defined in Glad and Ljung eq. 8.28. Note, this definition of closed loop is not the same as lft(P, K), which has B1 isntead of B2 as input matrix. Use `lft(l)` to get the system from disturbances to controlled variables `w -> z`.

The return value will be the closed loop from reference only, other disturbance signals (B1) are ignored. See `feedback` for a more advanced option.

Use `static_gain_compensation` to adjust the gain from references acting on the input B2, `dcgain(closedloop(l))*static_gain_compensation(l) ≈ I`
"""
function closedloop(l::LQGProblem, L = lqr(l), K = kalman(l))
    # todo: reimplement as lft
    P = system_mapping(l, identity)
    @unpack A, B1, B2, C2, C1 = l
    n = P.nx
    Acl = [A-B2*L B2*L; zero(A) A-K*C2] # 8.28
    #Glad Ljung has B2 here instead of B1. The difference lies in Glad, Ljung calling the system from references acting through B2 the closed loop, whereas most other literature uses lft(P, K) as the closed loop, i.e., from B1
    Bcl = [B2; zero(B2)]

    Ccl = [C1 zero(C1)]
    syscl = ss(Acl, Bcl, Ccl, 0, l.timeevol)
end

ControlSystems.lft(l::LQGProblem) = lft(l.sys, -observer_controller(l))



system_mapping(l::LQGProblem, args...) = system_mapping(l.sys, args...)

performance_mapping(l::LQGProblem, args...) = performance_mapping(l.sys, args...)


function Base.getproperty(G::LQGProblem, s::Symbol)
    if s ∈ fieldnames(LQGProblem)
        return getfield(G, s)
    end
    return getproperty(G.sys, s)
    error("No property named $s")
end

sensdoc = """
```
         ▲
         │e₁
         │  ┌─────┐
d₁────+──┴──►  P  ├─────┬──►e₄
      │     └─────┘     │
      │                 │
      │     ┌─────┐    -│
 e₂◄──┴─────┤  C  ◄──┬──+───d₂
            └─────┘  │
                     │e₃
                     ▼
```
- [`input_sensitivity`](@ref) is the transfer function from d₁ to e₁,       (I + CP)⁻¹
- [`output_sensitivity`](@ref) is the transfer function from d₂ to e₃,      (I + PC)⁻¹
- [`input_comp_sensitivity`](@ref) is the transfer function from d₁ to e₂,  (I + CP)⁻¹CP
- [`output_comp_sensitivity`](@ref) is the transfer function from d₂ to e₄, (I + PC)⁻¹PC
- [`G_PS`](@ref) is the transfer function from d₁ to e₄,                    (1 + PC)⁻¹P
- [`G_CS`](@ref) is the transfer function from d₂ to e₂,                    (1 + CP)⁻¹C
- [`feedback_control`](@ref) is the transfer function from d₂ to [e₄; e₂] (r → [y; u])
"""

"""
See [`output_sensitivity`](@ref)
$sensdoc
"""
function sensitivity(args...)# Sensitivity function
    return output_sensitivity(args...)
end

"""
See [`output_comp_sensitivity`](@ref)
$sensdoc
"""
function comp_sensitivity(args...) # Complementary sensitivity function
    return output_comp_sensitivity(args...)
end

"""
    G_PS(P, C)

The closed-loop transfer function from load disturbance to plant output.
Technically, it's `(1 + PC)⁻¹P` so `SP` would be a better, but nonstandard name.
$sensdoc
"""
G_PS(P, C) = output_sensitivity(P, C)*P
function G_PS(l::LQGProblem) # Load disturbance to output
    return output_sensitivity(l) * system_mapping(l)
end

"""
    G_CS(P, C)

The closed-loop transfer function from (-) measurement noise or (+) reference to control signal.
Technically, it's `(1 + CP)⁻¹C` so `SC` would be a better, but nonstandard name.
$sensdoc
"""
G_CS(P, C) = input_sensitivity(P, C)*C
function G_CS(l::LQGProblem) # Noise to control signal
    return input_sensitivity(l) * observer_controller(l)
end

# loopgain(P,C) = P*C
# function loopgain(l::LQGProblem)
#     return system_mapping(l)*observer_controller(l)
# end

# function returndifference(l::LQGProblem)
#     PC = loopgain(l)
#     p = size(l.C2, 1)
#     return ss(Matrix{numeric_type(PC)}(I, p, p), l.timeevol) + PC
# end

# function stabilityrobustness(l::LQGProblem)
#     PC = loopgain(l)
#     p = size(l.C2, 1)
#     return ss(Matrix{numeric_type(PC)}(I, p, p), l.timeevol) + inv(PC)
# end

"""
    input_sensitivity(P, C)
    input_sensitivity(l::LQGProblem)

Transfer function from load disturbance to total plant input.
- "Input" signifies that the transfer function is from the input of the plant.
$sensdoc
"""
function input_sensitivity(P,C)
    T = feedback(C * P)
    ss(I(noutputs(T)), P.timeevol) - T
end
input_sensitivity(l::LQGProblem) = input_sensitivity(system_mapping(l), observer_controller(l))

"""
    output_sensitivity(P, C)
    output_sensitivity(l::LQGProblem)

Transfer function from measurement noise / reference to control error.
- "output" signifies that the transfer function is from the output of the plant.
$sensdoc
"""
function output_sensitivity(P,C)
    PC = P*C
    S = feedback(ss(Matrix{numeric_type(PC)}(I, ninputs(PC), ninputs(PC)), P.timeevol), PC)
    S.C .*= -1
    S.B .*= -1
    S
end
output_sensitivity(l::LQGProblem) = output_sensitivity(system_mapping(l), observer_controller(l))

"""
    input_comp_sensitivity(P,C)
    input_comp_sensitivity(l::LQGProblem)

Transfer function from load disturbance to control signal.
- "Input" signifies that the transfer function is from the input of the plant.
- "Complimentary" signifies that the transfer function is to an output (in this case controller output)
$sensdoc
"""
function input_comp_sensitivity(P,C)
    T = feedback(C * P)
end
input_comp_sensitivity(l::LQGProblem) = input_comp_sensitivity(system_mapping(l), observer_controller(l))

"""
    output_comp_sensitivity(P,C)
    output_comp_sensitivity(l::LQGProblem)

Transfer function from measurement noise / reference to plant output.
- "output" signifies that the transfer function is from the output of the plant.
- "Complimentary" signifies that the transfer function is to an output (in this case plant output)
$sensdoc
"""
function output_comp_sensitivity(P,C)
    S = output_sensitivity(P,C)
    ss(I(noutputs(S)), P.timeevol) - S
end
output_comp_sensitivity(l::LQGProblem) = output_comp_sensitivity(system_mapping(l), observer_controller(l))

"""
    G = feedback_control(P, K)

Return the (negative feedback) closed-loop system from input of `K` to output of `P`
while outputing also the control signal (output of `K`), i.e.,
`G` maps references to `[y; u]`

# Example:
The following are two equivalent ways of achieving the same thing
```julia
G = ssrand(3,4,2)
K = ssrand(4,3,2)

Gcl1 = feedback_control(G, K) # First option

# Second option using named systems and connect
G = named_ss(G, :G)
K = named_ss(K, :K)
S = sumblock("Ku = r - Gy", n=3) # Create a sumblock that computes r - Gy for vectors of length 3

z1 = [G.y; K.y] # Output both plant and controller outputs
w1 = :r^3       # Extenal inputs are the three references into the sum block
connections = [K.y .=> G.u; G.y .=> G.y; K.u .=> K.u] # Since the sumblock uses the same names as the IO signals of G,K, we can reuse these names here
Gcl2 = connect([G, K, S], connections; z1, w1)

@test linfnorm(minreal(Gcl1 - Gcl2.sys))[1] < 1e-10 # They are the same
```

See also [`extended_gangoffour`](@ref).
"""
function feedback_control(G, K)
    ny,nu = size(G)
    Zperm = [(1:ny).+nu; 1:nu] # To make output come before control
    feedback(K, G; Z2 = :, Zperm)
end


plot(G::LQGProblem) = gangoffourplot(G)

function gangoffour(l::LQGProblem)
    sensitivity(l), G_PS(l), G_CS(l), comp_sensitivity(l)
end

function ControlSystems.gangoffourplot(l::LQGProblem, args...; sigma = true, kwargs...)
    S,D,N,T = gangoffour(l)
    bp = (args...; kwargs...) -> sigma ? sigmaplot(args...; kwargs...) : bodeplot(args...; plotphase=false, kwargs...)
    f1 = bp(S, args...; show=false, title="S = 1/(1+PC)", kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    try
        mag, freq = hinfnorm(S)
        isfinite(mag) && isfinite(freq) && (freq>0) && Plots.scatter!([freq], [mag], label="Mₛ = $(round(mag, digits=2))")
    catch
    end
    f2 = bodeplot(D, args...; show=false, title="D = P/(1+PC)", plotphase=false, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    f3 = bodeplot(N, args...; show=false, title="N = C/(1+PC)", plotphase=false, kwargs...)
    f4 = bp(T, args...; show=false, title="T = PC/(1+PC)", kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    try
        mag, freq = hinfnorm(T)
        isfinite(mag) && isfinite(freq) && (freq>0) && Plots.scatter!([freq], [mag], label="Mₜ = $(round(mag, digits=2))")
    catch
    end
    Plots.plot(f1,f2,f3,f4)
end


"""
    gangoffourplot2(P, C, args...; sigma = true, kwargs...)

Plot the gang-of-four. If the closed-loop is MIMO, output sensitivity functions will be shown.

$sensdoc
"""
function gangoffourplot2(P, C, args...; sigma = true, kwargs...)
    S,D,N,T = gangoffour2(P,C)
    bp = (args...; kwargs...) -> sigma ? sigmaplot(args...; kwargs...) : bodeplot(args...; plotphase=false, kwargs...)
    f1 = bp(S, args...; show=false, title="S = 1/(1+PC)", kwargs...)
    # Plots.hline!([1], l=(:black, :dash), primary=false)
    Plots.hline!([1.0 1.1 1.2], l=(:dash, [:green :orange :red]), sp=1, lab=["1.0" "1.1" "1.2"], ylims=(-3,1.8))
    f2 = bodeplot(D, args...; show=false, title="D = P/(1+PC)", plotphase=false, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    f3 = bodeplot(N, args...; show=false, title="N = C/(1+PC)", plotphase=false, kwargs...)
    f4 = bp(T, args...; show=false, title="T = PC/(1+PC)", ylims=(-3,1.8), kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    Plots.plot(f1,f2,f3,f4, ticks=:default, ylabel="", legend=:bottomright)
end

function gangofsevenplot(P, C, F, args...; sigma = true, ylabel="", layout=4, kwargs...)
    S,D,CS,T = gangoffour2(P,C)
    RY = T*F
    RU = CS*F
    RE = S*F
    bp = (args...; kwargs...) -> sigma ? sigmaplot!(args...; kwargs...) : bodeplot!(args...; plotphase=false, kwargs...)
    Plots.plot(; layout, ticks=:default, xscale=:log10, link=:both)
    bp(S, args...; show=false, title="S = 1/(1+PC)", sp=1, kwargs...)
    # Plots.hline!([1], l=(:black, :dash), primary=false, sp=1)
    Plots.hline!([1.0 1.1 1.2], l=(:dash, [:green :orange :red]), sp=1, lab=["1.0" "1.1" "1.2"], ylims=(-3,1.8))
    bodeplot!(D, args...; show=false, title="PS = P/(1+PC)", plotphase=false, sp=2, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false, sp=2)
    bodeplot!(CS, args...; show=false, title="CS = C/(1+PC)", plotphase=false, sp=3, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false, sp=3)
    bp(T, args...; show=false, title="T = PC/(1+PC)", sp=4, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false, sp=4)

    bp(RE, args...; show=false, title="S = 1/(1+PC)", lab="SF = r->e", l=(:dash, :red), sp=1, kwargs...)
    bp(RY, args...; show=false, title="T = PC/(1+PC)", lab="TF = r->y", l=(:dash, :green), sp=4, kwargs...)
    bp(RU, args...; show=false, title="CS = C/(1+PC)", lab="CSF = r->u", l=(:dash, :blue), sp=3, kwargs...)



    # Plots.plot(f1,f2,f3,f4,f7,f6,f5; ticks=:default, layout)
end

"""
    S, PS, CS, T = gangoffour2(P::LTISystem, C::LTISystem)

Return the gang-of-four. If the closed-loop is MIMO, output sensitivity functions will be shown.

$sensdoc
"""
function gangoffour2(P::LTISystem, C::LTISystem)
    S = output_sensitivity(P, C)
    PS = G_PS(P, C)
    CS = G_CS(P, C)
    T = output_comp_sensitivity(P, C)
    return S, PS, CS, T
end

