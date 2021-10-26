
import Base.getindex
import ControlSystems.numeric_type

"""
    G = LQG(sys::AbstractStateSpace, Q1, Q2, R1, R2; qQ=0, qR=0, M = I, N = I)

Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)`
where the controller is of LQG-type. The controller is specified by weight matrices `Q1,Q2`
that penalizes state deviations and control signal variance respectively, and covariance
matrices `R1,R2` which specify state drift and measurement covariance respectively.

`qQ` and `qR` can be set to incorporate loop transfer recovery, i.e.,
```julia
L = lqr(A, B, Q1+qQ*C'C, Q2)
K = kalman(A, C, R1+qR*B*B', R2)
```

`M` is a matrix that defines the controlled variables `z`, i.e., the variables for which you provide reference signals. If no `M` is provided, the default is to consider all state variables of the system as controlled. The definitions of `z` and `y` are given below
```
y = C*x
z = M*x
```
`size(M, 1)` determines the size of the `Q1` matrix you need to supply.

`N` is a matrix that defines how the dynamics noise `v` enters the system, i.e. If no `N` is provided, the default is to consider all state variables being affected by independent noise components. The definition of `v` is given below
```
x′ = A*x + B*u + N*v
```
`size(N, 2)` determines the size of the `R1` matrix you need to supply.

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

G = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR)

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
    kwargs...,
)
    sys isa ExtendedStateSpace || (sys = ExtendedStateSpace(sys))
    @unpack B1, C1, C2 = sys
    Q1 = Q1 isa AbstractVector ? diagm(Q1) : Q1
    Q2 = Q2 isa AbstractVector ? diagm(Q2) : Q2
    R1 = R1 isa AbstractVector ? diagm(R1) : R1
    R2 = R2 isa AbstractVector ? diagm(R2) : R2
    size(Q1, 1) == size(C1,1) || throw(ArgumentError("The size of Q1 is determined by C1, not by the state."))
    size(R2, 1) == size(C2,1) || throw(ArgumentError("The size of R2 is determined by C2, not by the state."))
    size(R1, 1) == size(B1,2) || throw(ArgumentError("The size of R1 is determined by B1, not by the state."))
    LQGProblem(sys, Q1, Q2, R1, R2, qQ, qR; kwargs...)
end 

"""
    LQGProblem(P::ExtendedStateSpace)

If only an `ExtendedStateSpace` system is provided, the system `P` is assumed to correspond to the H₂ optimal control problem with
```
C1'C1    = Q1
D12'D12  = Q2

B1*B1'   = R1
D21*D21' = R2
```
and an `LQGProblem` with the above covariance matrices is returned. The system
description in the returned LQGProblem will have `B1 = C1 = I`.
See Ch. 13 in Robust and optimal control for reference. 
"""
function LQGProblem(P::ExtendedStateSpace)
    @unpack A, B1, B2, C1, C2, D11, D12, D21, D22 = P
    # all(iszero, D22) || throw(ArgumentError("Non-zero D22 not handled."))
    all(iszero, D11) || throw(ArgumentError("Non-zero D11 not handled."))
    all(iszero, D12'C1) || throw(ArgumentError("D12'C1 should be 0"))
    all(iszero, D21*B1') || throw(ArgumentError("D21*B1' should be 0"))
    Q1 = C1'C1
    Q2 = D12'D12
    R1 = B1*B1'
    R2 = D21*D21'
    B1 = I(P.nx)
    C1 = I(P.nx)

    P = ss(A, B1, B2, C1, C2; D22, Ts = P.timeevol)
    # P = ss(A, B1, B2, C1, C2; Ts = P.timeevol)
    LQGProblem(P, Q1, Q2, R1, R2)
end

function ControlSystems.kalman(l::LQGProblem)
    @unpack A, C2, B1, R1, qR, B2, R2 = l
    fun = isdiscrete(l) ? dkalman : kalman
    K = fun(A, C2, Hermitian(B1*R1*B1' + qR * B2 * B2'), R2)
end

function ControlSystems.lqr(l::LQGProblem)
    @unpack A, B2, C1, Q1, qQ, C2, Q2 = l 
    fun = isdiscrete(l) ? dlqr : lqr
    L = fun(A, B2, Hermitian(C1'Q1*C1 + qQ * C2'C2), Q2)
end

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
```
C = observer_controller(l)
Ce = extended_controller(l)
system_mapping(Ce) == -C
````
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
    A,B,C,D = ssdata(system_mapping(l))
    Ac = A - B*L - K*C + K*D*L # 8.26b
    Bc = K
    Cc = L
    Dc = 0
    ss(Ac, Bc, Cc, Dc, l.timeevol)
end

function ff_controller(l::LQGProblem, L = lqr(l), K = kalman(l))
    Ae,Be,Ce,De = ssdata(system_mapping(l))
    Ac = Ae - Be*L - K*Ce + K*De*L # 8.26c
    Bc = Be*static_gain_compensation(l, L)
    Cc = L
    Dc = 0
    return 1 - ss(Ac, Bc, Cc, Dc, l.timeevol)
end

"""
    closedloop(l::LQGProblem, L = lqr(l), K = kalman(l))

Closed-loop system as defined in Glad and Ljung eq. 8.28

The return value will be the closed loop from reference only, other disturbance signals (B1) are ignored. See `feedback` for a more advanced option.
"""
function closedloop(l::LQGProblem, L = lqr(l), K = kalman(l))
    # todo: reimplement as lft
    P = system_mapping(l)
    @unpack A, B2, C2, C1 = l
    n = P.nx
    # Lr = pinv(C1 * ((P.B * L[:, 1:n] - P.A) \ P.B))
    Lr = static_gain_compensation(l, L[:, 1:n])
    # Lr = (D - (C - D*L) * inv(A - B*L) * B)
    if any(!isfinite, Lr) || all(iszero, Lr)
        @warn "Could not compensate for static gain automatically." Lr
        Lr = 1
    end
    Acl = [A-B2*L B2*L; zero(A) A-K*C2] # 8.28
    BLr = B2 * Lr # QUESTION: should be B1 here? Glad Ljung has B2
    Bcl = [BLr; zero(BLr)]
    Ccl = [C1 zero(C1)]
    syscl = ss(Acl, Bcl, Ccl, 0, l.timeevol)
end

# function closedloop(l::LQGProblem)
#     K = observer_controller(l)
#     Ke = extended_controller(K)
#     lft(l.sys, Ke)
# end

system_mapping(l::LQGProblem) = system_mapping(l.sys)

performance_mapping(l::LQGProblem) = performance_mapping(l.sys)


function Base.getproperty(G::LQGProblem, s::Symbol)
    if s ∈ fieldnames(LQGProblem)
        return getfield(G, s)
    end
    return getproperty(G.sys, s)
    error("No property named $s")
end

function sensitivity(args...)# Sensitivity function
    return output_sensitivity(args...)
end

function comp_sensitivity(args...) # Complementary sensitivity function
    return output_comp_sensitivity(args...)
end

G_PS(P, C) = P*input_sensitivity(P, C)
function G_PS(l::LQGProblem) # Load disturbance to output
    return system_mapping(l) * input_sensitivity(l)
end

G_CS(P, C) = C*output_sensitivity(P, C)
function G_CS(l::LQGProblem) # Noise to control signal
    return observer_controller(l) * output_sensitivity(l)
end

loopgain(P,C) = P*C
function loopgain(l::LQGProblem)
    return system_mapping(l)*observer_controller(l)
end

function returndifference(l::LQGProblem)
    PC = loopgain(L)
    p = size(l.C2, 1)
    return ss(Matrix{numeric_type(PC)}(I, p, p), l.timeevol) + PC
end

function stabilityrobustness(l::LQGProblem)
    PC = loopgain(L)
    p = size(l.C2, 1)
    return ss(Matrix{numeric_type(PC)}(I, p, p), l.timeevol) + inv(PC)
end

input_sensitivity(l::LQGProblem) = input_sensitivity(system_mapping(l), observer_controller(l))
function input_sensitivity(P,C)
    T = feedback(C * P)
    ss(I(noutputs(T)), P.timeevol) - T
end

input_comp_sensitivity(l::LQGProblem) = input_comp_sensitivity(system_mapping(l), observer_controller(l))
function input_comp_sensitivity(P,C)
    T = feedback(C * P)
end

output_sensitivity(l::LQGProblem) = output_sensitivity(system_mapping(l), observer_controller(l))
function output_sensitivity(P,C)
    PC = P*C
    S = feedback(ss(Matrix{numeric_type(PC)}(I, ninputs(PC), ninputs(PC)), P.timeevol), PC)
    S.C .*= -1
    S.B .*= -1
    S
end

output_comp_sensitivity(l::LQGProblem) = output_comp_sensitivity(system_mapping(l), observer_controller(l))
function output_comp_sensitivity(P,C)
    S = output_sensitivity(P,C)
    ss(I(noutputs(S)), P.timeevol) - S
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

function gangoffourplot2(P, C, args...; sigma = true, kwargs...)
    S,D,N,T = gangoffour2(P,C)
    bp = (args...; kwargs...) -> sigma ? sigmaplot(args...; kwargs...) : bodeplot(args...; plotphase=false, kwargs...)
    f1 = bp(S, args...; show=false, title="S = 1/(1+PC)", kwargs...)
    # Plots.hline!([1], l=(:black, :dash), primary=false)
    Plots.hline!([1.1 1.2 0.1], l=(:dash, [:orange :red :green]), sp=1, lab=["1.1" "1.2" "0.1"])
    f2 = bodeplot(D, args...; show=false, title="D = P/(1+PC)", plotphase=false, kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    f3 = bodeplot(N, args...; show=false, title="N = C/(1+PC)", plotphase=false, kwargs...)
    f4 = bp(T, args...; show=false, title="T = PC/(1+PC)", kwargs...)
    Plots.hline!([1], l=(:black, :dash), primary=false)
    Plots.plot(f1,f2,f3,f4, ticks=:default, ylabel="")
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
    Plots.hline!([1.1 1.2 0.1], l=(:dash, [:orange :red :green]), lab=["1.1" "1.2" "0.1"], sp=1)
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


function gangoffour2(P::LTISystem, C::LTISystem)
    S = output_sensitivity(P, C)
    PS = G_PS(P, C)
    CS = G_CS(P, C)
    T = output_comp_sensitivity(P, C)
    return S, PS, CS, T
end

