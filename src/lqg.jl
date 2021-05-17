
import Base.getindex
import ControlSystems.numeric_type

"""
    G = LQG(sys::AbstractStateSpace, Q1, Q2, R1, R2; qQ=0, qR=0, integrator=false, M = I, N = I)

Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)`
where the controller is of LQG-type. The controller is specified by weight matrices `Q1,Q2`
that penalizes state deviations and control signal variance respectively, and covariance
matrices `R1,R2` which specify state drift and measurement covariance respectively.
This constructor calls [`lqr`](@ref) and [`kalman`](@ref) and forms the closed-loop system.

If `integrator=true`, the resulting controller will have integral action.
This is achieved by adding a model of a constant disturbance on the inputs to the system
described by `sys`.

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

# Fields and properties
When the LQG-object is populated by the lqg-function, the following fields have been made available
- `L` is the feedback matrix, such that `A-BL` is stable. Note that the length of the state vector (and the width of L) is increased by the number of inputs if the option `integrator=true`.
- `K` is the kalman gain such that `A-KC` is stable

Several other properties of the object are accessible as properties. The available properties are
(some have many alternative names, separated with / )

- `G.cl / G.closedloop` is the closed-loop system, including observer, from reference to output, precompensated to have static gain 1 (`u = −Lx + lᵣr`).
- `G.S / G.Sin` Input sensitivity function
- `G.T / G.Tin` Input complementary sensitivity function
- `G.Sout` Output sensitivity function
- `G.Tout` Output complementary sensitivity function
- `G.CS` The transfer function from measurement noise to control signal
- `G.DS` The transfer function from input load disturbance to output
- `G.lt / G.looptransfer / G.loopgain  =  PC`
- `G.rd / G.returndifference  =  I + PC`
- `G.sr / G.stabilityrobustness  =  I + inv(PC)`
- `G.Fy / G.controller` Returns the controller as a StateSpace-system `u = L*inv(sI - A + BL + KC)*K * y`. This controller is acting on the measured signal, not the reference. The controller acting on the reference is `G.Fr`
- `G.Fr` Returns the controller from reference as a StateSpace-system. `I - L*inv(sI - A + BL + KC)*B`

It is also possible to access all fileds using the `G.symbol` syntax, the fields are `P,Q1,Q2,R1,R2,qQ,qR,sysc,L,K,integrator`

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

G = LQG(sys, Q1, Q2, R1, R2, qQ=qQ, qR=qR, integrator=true)

Gcl = G.cl
T = G.T
S = G.S
sigmaplot([S,T],exp10.(range(-3, stop=3, length=1000)))
stepplot(Gcl)
```

"""
struct LQGProblem
    P::ExtendedStateSpace
    Q1::AbstractMatrix
    Q2::AbstractMatrix
    R1::AbstractMatrix
    R2::AbstractMatrix
    qQ::Real
    qR::Real
end


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

function ControlSystems.kalman(l::LQGProblem)
    @unpack A, C2, B1, R1, qR, B2, R2 = l
    K = kalman(A, C2, B1*R1*B1' + qR * B2 * B2', R2)
end

function ControlSystems.lqr(l::LQGProblem)
    @unpack A, B2, C1, Q1, qQ, C2, Q2 = l 
    L = lqr(A, B2, C1'Q1*C1 + qQ * C2'C2, Q2)
end

function static_gain_compensation(l::LQGProblem, L = lqr(l))
    @unpack A, C1, B2 = l
    Lr = pinv(C1 * ((B2 * L - A) \ B2))
end


"""
    extended_controller(K)

Takes a controller and returns an `ExtendedStateSpace` version which has augmented input `[r; y]` and output `y` (`z` output is 0-dim).
"""
function extended_controller(K)
    nx,nu,ny = K.nx, K.nu, K.ny
    A,B,C,D = ssdata(K)
    ss(A, B, -B, zeros(0,nx), C, D21=D, D22=-D)
end

function controller(l::LQGProblem, L = lqr(l), K = kalman(l))
    A,B,C,D = ssdata(system_mapping(l))
    Ac = A - B*L - K*C + K*D*L # 8.26b
    Bc = K
    Cc = L
    Dc = 0
    ss(Ac, Bc, Cc, Dc)
end

function ff_controller(l::LQGProblem, L = lqr(l), K = kalman(l))
    Ae,Be,Ce,De = ssdata(system_mapping(l))
    Ac = Ae - Be*L - K*Ce + K*De*L # 8.26b
    Bc = Be*static_gain_compensation(l, L)
    Cc = L
    Dc = 0
    return 1 - ss(Ac, Bc, Cc, Dc)
end

function closedloop(l::LQGProblem, L = lqr(l), K = kalman(l))
    # todo: reimplement as lft
    P = system_mapping(l)
    @unpack A, B2, C2, C1 = l
    n = P.nx
    Lr = pinv(C1 * ((P.B * L[:, 1:n] - P.A) \ P.B))
    if any(!isfinite, Lr) || all(iszero, Lr)
        @warn "Could not compensate for static gain automatically." Lr
        Lr = 1
    end
    Acl = [A-B2*L B2*L; zero(A) A-K*C2] # 8.28
    BLr = B2 * Lr
    Bcl = [BLr; zero(BLr)]
    Ccl = [C1 zero(C1)]
    syscl = ss(Acl, Bcl, Ccl, 0)
end

# function closedloop(l::LQGProblem)
#     K = controller(l)
#     Ke = extended_controller(K)
#     lft(l.P, Ke)
# end

system_mapping(l::LQGProblem) = system_mapping(l.P)

performance_mapping(l::LQGProblem) = performance_mapping(l.P)


function Base.getproperty(G::LQGProblem, s::Symbol)
    if s ∈ fieldnames(LQGProblem)
        return getfield(G, s)
    end
    if s ∈ (:A, :B, :C, :D, :timeevol, :Ts, :B1, :B2, :C1, :C2, :D11, :D12, :D21, :D22)
        return getproperty(G.P, s)
    end
    error("No property named $s")
end

function sensitivity(l::LQGProblem)# Sensitivity function
    return output_sensitivity(l)
end

function comp_sensitivity(l::LQGProblem) # Complementary sensitivity function
    return output_comp_sensitivity(l)
end

function G_PS(l::LQGProblem) # Load disturbance to output
    return system_mapping(l) * input_sensitivity(l)
end

function G_CS(l::LQGProblem) # Noise to control signal
    return controller(l) * output_sensitivity(l)
end

function loopgain(l::LQGProblem)
    return system_mapping(l)*controller(l)
end

function returndifference(l::LQGProblem)
    PC = loopgain(L)
    p = size(l.C2, 1)
    return ss(Matrix{numeric_type(PC)}(I, p, p)) + PC
end

function stabilityrobustness(l::LQGProblem)
    PC = loopgain(L)
    p = size(l.C2, 1)
    return ss(Matrix{numeric_type(PC)}(I, p, p)) + inv(PC)
end

input_sensitivity(l::LQGProblem) = input_sensitivity(system_mapping(l), controller(l))
function input_sensitivity(P,C)
    T = feedback(C * P)
    ss(I(noutputs(T))) - T
end

input_comp_sensitivity(l::LQGProblem) = input_comp_sensitivity(system_mapping(l), controller(l))
function input_comp_sensitivity(P,C)
    T = feedback(C * P)
end

output_sensitivity(l::LQGProblem) = output_sensitivity(system_mapping(l), controller(l))
function output_sensitivity(P,C)
    PC = P*C
    S = feedback(ss(Matrix{numeric_type(PC)}(I, ninputs(PC), ninputs(PC))), PC)
    S.C .*= -1
    S.B .*= -1
    S
end

output_comp_sensitivity(l::LQGProblem) = output_comp_sensitivity(system_mapping(l), controller(l))
function output_comp_sensitivity(P,C)
    S = output_sensitivity(P,C)
    ss(I(noutputs(S))) - S
end


plot(G::LQGProblem) = gangoffourplot(G)

function gangoffour(l::LQGProblem)
    sensitivity(l), G_PS(l), G_CS(l), comp_sensitivity(l)
end

function ControlSystems.gangoffourplot(l::LQGProblem, args...; kwargs...)
    S,D,N,T = gangoffour(l)
    f1 = sigmaplot(S, args...; show=false, title="\$S = 1/(1+PC)\$", kwargs...)
    f2 = sigmaplot(D, args...; show=false, title="\$D = P/(1+PC)\$", kwargs...)
    f3 = sigmaplot(N, args...; show=false, title="\$N = C/(1+PC)\$", kwargs...)
    f4 = sigmaplot(T, args...; show=false, title="\$T = PC/(1+PC\$)", kwargs...)
    Plots.plot(f1,f2,f3,f4)
end


