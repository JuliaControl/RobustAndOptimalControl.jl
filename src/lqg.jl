
import Base.getindex
import ControlSystemsBase.numeric_type

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
In this example we will control a MIMO system with one unstable pole and one unstable zero. When the system has both unstable zeros and poles, there are fundamental limitations on performance. The unstable zero is in this case faster than the unstable pole, so the system is controllable. For good performance, we want as large separation between the unstable zero dynamics and the unstable poles as possible. 
```julia
s = tf("s")
P = [1/(s+1) 2/(s+2); 1/(s+3) 1/(s-1)]
sys = ExtendedStateSpace(ss(P)) # Controlled outputs same as measured outputs and state noise affects at inputs only. 
eye(n) = Matrix{Float64}(I,n,n) # For convinience

qQ = 0
qR = 0
Q1 = 10eye(2)
Q2 = 1eye(2)
R1 = 1eye(2)
R2 = 1eye(2)

G = LQGProblem(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR)

T = comp_sensitivity(G)
S = sensitivity(G)
Gcl = closedloop(G)*static_gain_compensation(G)
plot(
    sigmaplot([S,T, Gcl],exp10.(range(-3, stop=3, length=1000)), lab=["S" "T" "Gry"]),
    plot(step(Gcl, 5))
)
```

# Extended help
Several functions are defined for instances of `LQGProblem`
- [`closedloop`](@ref)
- [`extended_controller`](@ref)
- [`ff_controller`](@ref)
- [`gangoffour`](@ref)
- [`G_CS`](@ref)
- [`G_PS`](@ref)
- [`input_comp_sensitivity`](@ref)
- [`input_sensitivity`](@ref)
- [`output_comp_sensitivity`](@ref)
- [`output_sensitivity`](@ref)
- [`system_mapping`](@ref)
- [`performance_mapping`](@ref)
- [`static_gain_compensation`](@ref)
- [`gangoffourplot`](@ref)
- [`kalman`](@ref)
- [`lft`](@ref)
- [`lqr`](@ref)
- [`observer_controller`](@ref)

A video tutorial on how to use the LQG interface is available [here](https://youtu.be/NuAxN1mGCPs)
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

ControlSystemsBase.isdiscrete(l::LQGProblem) = ControlSystemsBase.isdiscrete(l.sys)
ControlSystemsBase.iscontinuous(l::LQGProblem) = ControlSystemsBase.iscontinuous(l.sys)

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
    if sys isa ExtendedStateSpace
        iszero(sys.D12) || error("When providing quadratic penatly matrices, non-zero D12 is not supported since it is not clear how to interpret the control input cost. Consider using the H2 interface instead.")
        iszero(sys.D21) || error("When providing covariance matrices, non-zero D21 is not supported since it is not clear how to interpret the measurement noise. Consider using the H2 interface instead.")
    else
        sys = ExtendedStateSpace(sys, B1=I, C1=I)
    end
    @unpack B1, C1, C2 = sys
    Q1 = Q1 isa AbstractVector ? diagm(Q1) : Q1
    Q2 = Q2 isa AbstractVector ? diagm(Q2) : Q2
    R1 = R1 isa AbstractVector ? diagm(R1) : R1
    R2 = R2 isa AbstractVector ? diagm(R2) : R2
    size(Q1, 1) == size(C1,1) || throw(ArgumentError("The size of Q1 is determined by C1, not by the state, expected size = $(size(C1,1))."))
    size(R2, 1) == size(C2,1) || throw(ArgumentError("The size of R2 is determined by C2, not by the state, expected size = $(size(C2,1))."))
    size(R1, 1) == size(B1,2) || throw(ArgumentError("The size of R1 is determined by B1, not by the state, expected size = $(size(B1,2))."))
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

function ControlSystemsBase.kalman(l::LQGProblem; direct = false)
    @unpack A, C2, B1, R1, qR, B2, R2, SR = l
    # We do not apply the transformation D21*R2*D21' since when the user has provided an ESS, R2 == D21*D21', and when the user provides covariance matrices, R2 is provided directly.
    K = kalman(l.timeevol, A, C2, Hermitian(B1*R1*B1' + qR * B2 * B2'), R2, SR; direct)
end

function ControlSystemsBase.lqr(l::LQGProblem)
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
function lqr3(P::AbstractStateSpace{<:Discrete}, Q1::AbstractMatrix, Q2::AbstractMatrix, Q3::AbstractMatrix; full=false)
    Pd = add_input_differentiator(P)
    S = zeros(Pd.nx, P.nu)
    S[P.nx+1:end, :] = -Q3
	X, _, L = MatrixEquations.ared(Pd.A, Pd.B, Q2+Q3, cat(Q1, Q3, dims=(1,2)), S) # ME has cost matrices reversed
    full ? L : L[:, 1:P.nx]
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
    # The Δu term expands to u⁺'Q3u⁺ + u'Q3u - 2u Q3 u⁺, so the factor 2 is already accounted for
    Pd = add_input_differentiator(P)
    S = zeros(Pd.nx, P.nu)
    S[P.nx+1:end, :] = -Q3
	X, _, L = MatrixEquations.ared(Pd.A, Pd.B, Q2+Q3, cat(Q1, Q3, dims=(1,2)), S) # ME has cost matrices reversed
    full ? X : X[1:P.nx, 1:P.nx]
end

dare3(A::AbstractMatrix, B, Q1, Q2, Q3::AbstractMatrix) = dare3(ss(A, B, I(size(A,1)), 0, 1), Q1, Q2, Q3)



"""
    static_gain_compensation(l::LQGProblem, L = lqr(l))
    static_gain_compensation(A, B, C, D, L)

Find ``L_r`` such that
```
dcgain(closedloop(G)*Lr) ≈ I
```
"""
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

Returns an expression for the controller that is obtained when state-feedback `u = L(xᵣ-x̂)` is combined with a Kalman filter with gain `K` that produces state estimates x̂. The controller is an instance of `ExtendedStateSpace` where `C2 = -L, D21 = L` and `B2 = K`.

The returned system has *inputs* `[xᵣ; y]` and outputs the control signal `u`. If a reference model `R` is used to generate state references `xᵣ`, the controller from `(ry, y) -> u`  where `ry - y = e` is given by
```julia
Ce = extended_controller(l)
Ce = named_ss(ss(Ce); x = :xC, y = :u, u = [R.y; :y^l.ny]) # Name the inputs of Ce the same as the outputs of `R`.
connect([R, Ce]; u1 = R.y, y1 = R.y, w1 = [R.u; :y^l.ny], z1=[:u])
```

Since the negative part of the feedback is built into the returned system, we have
```julia
C = observer_controller(l)
Ce = extended_controller(l)
system_mapping(Ce) == -C
```

Please note, without the reference pre-filter, the DC gain from references to controlled outputs may not be identity.
"""
function extended_controller(l::LQGProblem, L = lqr(l), K = kalman(l))
    P = system_mapping(l)
    A,B,C,D = ssdata(P)
    Ac = A - B*L - K*C + K*D*L # 8.26b
    nx = l.nx
    B1 = zeros(nx, nx) # dynamics not affected by r
    # l.D21 does not appear here, see comment in kalman
    B2 = K # input y
    D21 = L #   L*xᵣ # should be D21?
    C2 = -L # - L*x̂
    C1 = zeros(0, nx)
    ss(Ac, B1, B2, C1, C2; D21, Ts = l.timeevol)
end

"""
    observer_controller(l::LQGProblem, L = lqr(l), K = kalman(l))

Returns an expression for the feedback controller ``C_{fb}`` in
``u = C_{fb}y + C_{ff}r``
that is obtained when state-feedback `u = -Lx̂` is combined with a Kalman filter with gain `K` that produces state estimates x̂.

Note: the transfer function returned is only a representation of the controller in the simple setting described above, e.g., it is not valid if the actual input contains anything that is not produced by a pure feedback from observed states. To obtain a controller that takes references into account, see `extended_controller`.

See also [`ff_controller`](@ref) that generates ``C_{ff}``.
"""
function ControlSystemsBase.observer_controller(l::LQGProblem, L::AbstractMatrix = lqr(l), K::Union{AbstractMatrix, Nothing} = nothing; direct = false)
    A,B,C,D = ssdata(system_mapping(l, identity))
    if K === nothing
        K = kalman(l; direct)
    end
    if direct && isdiscrete(sys)
        iszero(D) || throw(ArgumentError("D must be zero when using direct formulation of `observer_controller`"))
        IKC = (I - K*C)
        ABL = (A - B*L)
        Ac = IKC*ABL
        Bc = IKC*ABL*K
        Cc = L
        Dc = L*K
    else
        Ac = A - B*L - K*C + K*D*L # 8.26b
        Bc = K
        Cc = L
        Dc = 0
        iszero(l.D11) || error("Nonzero D11 not supported")
        iszero(l.D22) || error("Nonzero D22 not supported. The _transformP2Pbar is not used for LQG, but perhaps shpuld be?")
    end
    # do we need some way to specify which non-controllable inputs are measurable? No, because they will automatically appear in the measured outputs :)
    ss(Ac, Bc, Cc, Dc, l.timeevol)
end

"""
    ff_controller(l::LQGProblem, L = lqr(l), K = kalman(l))

Return the feedforward controller ``C_{ff}`` that maps references to plant inputs:
``u = C_{fb}y + C_{ff}r``

See also [`observer_controller`](@ref).
"""
function ff_controller(l::LQGProblem, L = lqr(l), K = kalman(l); comp_dc = true)
    Ae,Be,Ce,De = ssdata(system_mapping(l, identity))
    Ac = Ae - Be*L - K*Ce + K*De*L # 8.26c
    Cc = L
    Dc = 0
    if comp_dc
        Lr = static_gain_compensation(l, L)
        Bc = Be*Lr
        return Lr - ss(Ac, Bc, Cc, Dc, l.timeevol)
    else
        Bc = Be
        return I(size(Cc, 1)) - ss(Ac, Bc, Cc, Dc, l.timeevol)
    end
end
end

"""
    closedloop(l::LQGProblem, L = lqr(l), K = kalman(l))

Closed-loop system as defined in Glad and Ljung eq. 8.28. Note, this definition of closed loop is not the same as lft(P, K), which has B1 instead of B2 as input matrix. Use `lft(l)` to get the system from disturbances to controlled variables `w -> z`.

The return value will be the closed loop from filtred reference only, other disturbance signals (B1) are ignored. See [`feedback`](@ref) for a more advanced option. This function assumes that the control signal is computed as `u = r̃ - Lx̂` (not `u = L(xᵣ - x̂)`), i.e., the feedforward signal `r̃` is added directly to the plant input. `r̃` must thus be produced by an inverse-like model that takes state references and output the feedforward signal.

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

ControlSystemsBase.lft(l::LQGProblem) = lft(l.sys, -observer_controller(l))



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
    G_PS(l::LQGProblem)
"""
function G_PS(l::LQGProblem) # Load disturbance to output
    return output_sensitivity(l) * system_mapping(l)
end

"""
    G_CS(l::LQGProblem)
"""
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
    input_sensitivity(l::LQGProblem)
"""
input_sensitivity(l::LQGProblem) = input_sensitivity(system_mapping(l), observer_controller(l))

"""
    output_sensitivity(l::LQGProblem)
"""
output_sensitivity(l::LQGProblem) = output_sensitivity(system_mapping(l), observer_controller(l))

"""
    input_comp_sensitivity(l::LQGProblem)
"""
input_comp_sensitivity(l::LQGProblem) = input_comp_sensitivity(system_mapping(l), observer_controller(l))

"""
    output_comp_sensitivity(l::LQGProblem)
"""
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

To include also an input disturbance, use
```
Gcl = feedback(K, P, W2=:, Z2=:, Zperm=[(1:ny).+nu; 1:nu]) # y,u from r,d
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

function ControlSystemsBase.gangoffourplot(l::LQGProblem, args...; sigma = true, kwargs...)
    plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
    haskey(Base.loaded_modules, plots_id) || error("Call using Plots before calling this function")
    Plots = Base.loaded_modules[plots_id]
    S,D,N,T = gangoffour(l)
    bp = (args...; kwargs...) -> sigma ? sigmaplot(args...; kwargs...) : bodeplot(args...; plotphase=false, kwargs...)
    f1 = bp(S, args...; show=false, title="S", kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    try
        mag, freq = hinfnorm(S)
        isfinite(mag) && isfinite(freq) && (freq>0) && Plots.scatter!([freq], [mag], label="Mₛ = $(round(mag, digits=2))")
    catch
    end
    f2 = bodeplot(D, args...; show=false, plot_title="PS", plotphase=false, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    f3 = bodeplot(N, args...; show=false, plot_title="CS", plotphase=false, kwargs...)
    f4 = bp(T, args...; show=false, title="T", kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    try
        mag, freq = hinfnorm(T)
        isfinite(mag) && isfinite(freq) && (freq>0) && Plots.scatter!([freq], [mag], label="Mₜ = $(round(mag, digits=2))")
    catch
    end
    Plots.plot(f1,f2,f3,f4)
end

function gangofsevenplot(P, C, F, args...; sigma = true, ylabel="", layout=4, kwargs...)
    plots_id = Base.PkgId(UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
    haskey(Base.loaded_modules, plots_id) || error("Call using Plots before calling this function")
    Plots = Base.loaded_modules[plots_id]
    S,D,CS,T = gangoffour(P,C)
    RY = T*F
    RU = CS*F
    RE = S*F
    bp = (args...; kwargs...) -> sigma ? sigmaplot!(args...; kwargs...) : bodeplot!(args...; plotphase=false, kwargs...)
    Plots.plot(; layout, ticks=:default, xscale=:log10, link=:both)
    bp(S, args...; show=false, title="S = 1/(1+PC)", sp=1, kwargs...)
    # Plots.hline!([1], l=(:black, :dash), primary=false, sp=1)
    Plots.hline!([1.0 1.1 1.2], l=(:dash, [:green :orange :red]), sp=1, lab=["1.0" "1.1" "1.2"], ylims=(1e-3,8e1))
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

Base.@deprecate gangoffour2(P, C) gangoffour(P, C)
Base.@deprecate gangoffourplot2(P, C) gangoffourplot(P, C)

