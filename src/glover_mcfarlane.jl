@doc raw"""
    K, γ, info = glover_mcfarlane(G::AbstractStateSpace, γ = 1.1; W1=1, W2=1)

Design a controller for `G` that maximizes the stability margin ϵ = 1/γ with normalized coprime factor uncertainty using the method of Glover and McFarlane
```
γ = 1/ϵ = ||[K;I] inv(I-G*K)*inv(M)||∞
G = inv(M + ΔM)*(N + ΔN)
```
γ is given as a relative factor above γmin and must be greater than 1, i.e., if γ = 1.1, the controller will be designed for γ = 1.1*γmin.

We want γmin (which is always ≥ 1) as small as possible, and we usually require that γmin is less than 4, corresponding to 25% allowed coprime uncertainty.

Performance modeling is incorporated in the design by calling `glover_mcfarlane` on the shaped system `Gs = W2*G*W1` and then forming the controller as `K = W1*Ks*W2`. Using this formulation, traditional loop shaping can be done on `Gs = W2*G*W1`. The plant shaping is handled internally if keyword arguments `W1, W2` are used and the returned controller is already scaled. In this case, `Gs` and `Ks` are included in the `info` named tuple for inspection. 

See also [`glover_mcfarlane_2dof`](@ref) to design a feedforward filter as well.

# Example:
Example 9.3 from the reference below.
```julia
using RobustAndOptimalControl, ControlSystems, Plots, Test
G = tf(200, [10, 1])*tf(1, [0.05, 1])^2     |> ss
Gd = tf(100, [10, 1])                       |> ss
W1 = tf([1, 2], [1, 1e-6])                  |> ss
K, γ, info = glover_mcfarlane(G, 1.1; W1)
@test info.γmin ≈ 2.34 atol=0.005
Gcl = extended_gangoffour(G, K) # Form closed-loop system

bodeplot([G, info.Gs, G*K], lab=["G" "" "G scaled" "" "Loop transfer"]) |> display
bodeplot(Gcl, lab=["S" "KS" "PS" "T"], plotphase=false) |> display # Plot gang of four

plot( step(Gd*feedback(1, info.Gs), 3), lab="Initial controller")
plot!(step(Gd*feedback(1, G*K), 3), lab="Robustified") |> display

nyquistplot([info.Gs, G*K], ylims=(-2,1), xlims=(-2, 1),
    Ms_circles = 1.5,
    lab = ["Initial controller" "Robustified"],
    title = "Loop transfers with and without robustified controller"
    ) |> display
```

Ref: Sec 9.4.1 of Skogestad, "Multivariable Feedback Control: Analysis and Design"

# Extended help

Skogestad gives the following general advice:
1. Scale the plant outputs and inputs. This is very important for most design
    procedures and is sometimes forgotten. In general, scaling improves the
    conditioning of the design problem, it enables meaningful analysis to be made
    of the robustness properties of the feedback system in the frequency domain,
    and for loop-shaping it can simplify the selection of weights. There are a variety
    of methods available including normalization with respect to the magnitude of
    the maximum or average value of the signal in question. If one is to go straight to a design the following variation has
    proved useful in practice:
    (a) The outputs are scaled such that equal magnitudes of cross-coupling into each
        of the outputs is equally undesirable.
    (b) Each input is scaled by a given percentage (say 10%) of its expected range
        of operation. That is, the inputs are scaled to reflect the relative actuator
        capabilities. An example of this type of scaling is given in the aero-engine
        case study of Chapter 12.
2. Order the inputs and outputs so that the plant is as diagonal as possible. The
    relative gain array [`rga`](@ref) can be useful here. The purpose of this pseudo-diagonalization
    is to ease the design of the pre- and post-compensators which, for simplicity, will
    be chosen to be diagonal.

    Next, we discuss the selection of weights to obtain the shaped plant $G_s = W_2 G W_1$
    where $W_1 = W_p W_a W_g$
3. Select the elements of diagonal pre- and post-compensators $W_p$ and $W_2$ so that
    the singular values of $W_2 G W_p$ are desirable. This would normally mean high
    gain at low frequencies, roll-off rates of approximately 20 dB/decade (a slope of
    about 1) at the desired bandwidth(s), with higher rates at high frequencies. Some
    trial and error is involved here. $W_2$ is usually chosen as a constant, reflecting the
    relative importance of the outputs to be controlled and the other measurements
    being fed back to the controller. For example, if there are feedback measurements
    of two outputs to be controlled and a velocity signal, then $W_2$ might be chosen
    to be `diag([1, 1, 0.1])`, where 0.1 is in the velocity signal channel. $W_p$ contains the
    dynamic shaping. Integral action, for low frequency performance; phase-advance
    for reducing the roll-off rates at crossover, and phase-lag to increase the roll-off
    rates at high frequencies should all be placed in $W_p$ if desired. The weights should
    be chosen so that no unstable hidden modes are created in $G_s$.
5. Optional: Introduce an additional gain matrix $W_g$ cascaded with $W_a$ to provide
    control over actuator usage. $W_g$ is diagonal and is adjusted so that actuator rate
    limits are not exceeded for reference demands and typical disturbances on the
    scaled plant outputs. This requires some trial and error.
    
6. Robustly stabilize the shaped plant $G_s = W_2 G W_1$ , where $W_1 = W_p W_a W_g$,
    using `glover_mcfarlane`. First, the maximum stability
    margin $ϵ_{max} = 1/γ_{min}$ is calculated. If the margin is too small, $ϵmax < 0.25$, then go back and modify the weights. Otherwise, a γ-suboptimal controller is synthesized. There is usually no advantage to be gained by using the optimal controller. When $ϵ_{max}$ > 0.25
    (respectively $γ_{min}$ < 4) the design is usually successful. In this case, at least
    25% coprime factor uncertainty is allowed, and we also find that the shape of the
    open-loop singular values will not have changed much after robust stabilization.
    A small value of ϵmax indicates that the chosen singular value loop-shapes are
    incompatible with robust stability requirements. That the loop-shapes do not
    change much following robust stabilization if γ is small (ϵ large), is justified
    theoretically in McFarlane and Glover (1990).

7. Analyze the design and if all the specifications are not met make further
    modifications to the weights.
8. Implement the controller. The configuration shown in below has been found
    useful when compared with the conventional set up. This is because
    the references do not directly excite the dynamics of $K$, which can result in large amounts of overshoot (classical derivative kick). The constant prefilter ensures a steady-state gain of 1 between r and y, assuming integral action in $W_1$ or $G$ (note, the K returned by this function has opposite sign compared to that of Skogestad, so we use negative feedback here).

Anti-windup can be added to $W_1$ but putting $W_1$ on Hanus form after the synthesis, see [`hanus`](@ref).

```
       ┌─────────┐      ┌────────┐      ┌────────┐
    r  │         │    us│        │  u   │        │  y
   ───►│(K*W2)(0)├──+──►│   W1   ├─────►│   G    ├────┬──►
       │         │  │-  │        │      │        │    │
       └─────────┘  │   └────────┘      └────────┘    │
                    │                                 │
                    │                                 │
                    │   ┌────────┐      ┌────────┐    │
                    │   │        │  ys  │        │    │
                    └───┤   K    │◄─────┤   W2   │◄───┘
                        │        │      │        │
                        └────────┘      └────────┘
```
Keywords: `nfcsyn`, `coprimeunc`
"""
function glover_mcfarlane(G::AbstractStateSpace{Continuous}, γ = 1.1; W1=1, W2=1)
    γ > 1 || throw(ArgumentError("γ must be greater than 1"))
    Gs = W2*G*W1
    A,B,C,D = ssdata(Gs)

    R = I + D*D'
    S = I + D'D
    Sl = lu(S)
    # arec(A, B, R, Q, S) solves A'X + XA - (XB+S)R^(-1)(B'X+S') + Q = 0
    
    Ā = A - B*(S\D'C)
    Z,_ = MatrixEquations.arec(Ā', C', R, B*(Sl\B'))
    X,_ = MatrixEquations.arec(Ā, B, S, C'*(R\C))

    γmin = sqrt(1 + ρ(X*Z))

    γ *= γmin

    L = (1-γ^2)*I + X*Z
    F = -(Sl\(D'C + B'X))
    BK = γ^2*(L'\Z)*C'
    AK = A + B*F + BK*(C + D*F)
    CK = B'X
    DK = -D'
    Ks = -ss(AK, BK, CK, DK)
    Gcl = extended_gangoffour(Gs, Ks)
    imargin, ω = hinfnorm2(Gcl)
    K = W1*Ks*W2
    Gcl = extended_gangoffour(G, K)
    K, γ, (; Gcl, margin = inv(imargin), ω, γmin, Ks, Gs, Z, X)
end

"""
    K, γ, info = glover_mcfarlane(G::AbstractStateSpace{<:Discrete}, γ = 1.1; W1=1, W2=1, strictly_proper=false)

For discrete systems, the `info` tuple contains also feedback gains `F, L` and observer gain `Hkf` such that the controller on observer form is given by
```math
x^+ = Ax + Bu + H_{kf}*(Cx - y)\\\\
u = Fx + L*(Cx - y)
```
Note, this controller is *not* strictly proper, i.e., it has a non-zero D matrix.
The controller can be transformed to observer form for the scaled plant (`info.Gs`)
by `Ko = observer_controller(info)`, in which case the following holds `G*K == info.Gs*Ko` (realizations are different).

If `strictly_proper = true`, the returned controller `K` will have `D == 0`.
This can be advantageous in implementations where computational delays are present. In this case, `info.L == 0` as well.

Ref discrete version: Iglesias, "The Strictly Proper Discrete-Time Controller for the Normalized Left-Coprime Factorization Robust Stabilization Problem"
"""
function glover_mcfarlane(G::AbstractStateSpace{<:Discrete}, γ = 1.1; W1=1, W2=1, strictly_proper=false)
    γ > 1 || throw(ArgumentError("γ must be greater than 1"))
    Gs = W2*G*W1
    A,B,C,D = ssdata(Gs)
    iszero(D) || throw(ArgumentError("System must be strictly proper (D must be 0)"))

    X,_,Flq = MatrixEquations.ared(A, B, I, C'C)
    Z,_ = MatrixEquations.ared(A', C', I, B*B')
    Q1 = I + C*Z*C'

    if strictly_proper
        T = Q1 # They appear to be the same
        T12 = sqrt(T)
        X12 = sqrt(X)
        M1 = T12\C*Z*A'X12
        Q2 = I + X12*Z*X12

        QM = [
            sqrt(Q1 + 1/4 * M1*M1') -1/2*M1
            -1/2*M1'       sqrt(Q2 + 1/4 * M1'*M1)
        ]
        γmin = ρ(QM)
        γ *= γmin


        W = (γ^2 - 1)*I - Z*X
        XW = X/W
        Hkf = -A*Z*C'/Q1

        # Flq is the feedback gain from the ARE above, ref Rowe and Maciejowski, "Tuning MPC using H∞ Loop Shaping" between eq 36-37.
        # Flq = -B'X*((I + B*B'X)\A) # NOTE: Flq might be required for inverse optimal control
        F = -γ^2*Flq/W # Original paper
        # F = -γ^2*B'XW*((I + γ^2*B*B'XW)\A) # Merl paper
        # @show [F; Flq; F0]
        Ak = A + B*F + Hkf*C
        # NOTE: It's unclear what the best realization is. Both papers talk about a "dual" realization, but they appear to disagree on the formulation, and none of the expressions for the Q-matrices, in the Merl paper or the Maciejowski paper, lead to the same caluclated feedback gain. I think they might have left out a similarity transform on some matrices. The equations used here are thus not any of the dual forms, even if those appear to be recommended.

        L = 0
        Ck = F
        Bk = -Hkf
        Dk = 0
    else

        γmin = sqrt(1 + ρ(X*Z))
        γ *= γmin

        W = (γ^2 - 1)*I - Z*X
        Hkf = -A*Z*C'/Q1
        Akf = A+Hkf*C

        XW = X/W
        hest = lu(I + γ^2*B*B'XW)
        γ2B = γ^2*B'
        arne = γ2B*XW
        Ak = hest\Akf
        Ck = arne*Ak
        Bk = hest\Hkf
        Dk = arne*Bk

        F0 = -γ2B*XW/hest
        F = F0*A
        L = F0*Hkf
    end

    
    Ks = -ss(Ak, Bk, Ck, Dk, G.timeevol)
    Gcl = extended_gangoffour(Gs, Ks)
    imargin, ω = hinfnorm2(Gcl)
    K = W1*Ks*W2
    Gcl = extended_gangoffour(G, K)
    K, γ, (; Gcl, margin = inv(imargin), ω, γmin, Ks, Gs, Z, X, F, L, Hkf, W, W1, W2)
end

"""
    observer_controller(glover_mcfarlane_info::NamedTuple)

Return a controller on observer form (observer of the scaled plant `info.Gs`).
The observer controller has input vector `y`.

Ref: eq (2.5)-(2.6) of Iglesias, "The Strictly Proper Discrete-Time Controller for the Normalized Left-Coprime Factorization Robust Stabilization Problem"
"""
function ControlSystems.observer_controller(info::NamedTuple)
    isdiscrete(info.Gs) || throw(ArgumentError("Observer controller can only be generated for Glover McFarlane designs on discrete systems."))
    A,B,C,D = ssdata(info.Gs)
    iszero(D) || throw(ArgumentError("observer_controller does not support non-zero D matrix")) # not sure if this is a strict limitation or the paper author simplified.
    H, L, F = info.Hkf, info.L, info.F
    Ao = A + B*(F + L*C) + H*C
    Bo = (B*L+H)
    Co = F + L*C
    Do = L
    ss(Ao, Bo, Co, Do, info.Gs.timeevol)
end

"""
    observer_predictor(glover_mcfarlane_info::NamedTuple)

Return a predictor for the scaled plant `info.Gs`.
The observer predictor has input vector `[u; y]`.

Ref: eq (2.5)-(2.6) of Iglesias, "The Strictly Proper Discrete-Time Controller for the Normalized Left-Coprime Factorization Robust Stabilization Problem"
"""
function ControlSystems.observer_predictor(info::NamedTuple)
    isdiscrete(info.Gs) || throw(ArgumentError("Observer predictor can only be generated for Glover McFarlane designs on discrete systems."))
    A,B,C,D = ssdata(info.Gs)
    H, L, F = info.Hkf, info.L, info.F
    Ao = A + H*C
    Bo = [B-H*D -H]
    Co = C
    Do = [D zeros(size(D,1), size(H, 2))]
    ss(Ao, Bo, Co, Do, info.Gs.timeevol)
end

"""
    K, γ, info = glover_mcfarlane_2dof(G::AbstractStateSpace{Continuous}, Tref::AbstractStateSpace{Continuous}, γ = 1.1, ρ = 1.1;
    W1 = 1, Wo = I, match_dc = true, kwargs...)

Joint design of feedback and feedforward compensators
```math
K = \\left[K_1 & K_2\\right]
```
```
   ┌──────┐   ┌──────┐        ┌──────┐    ┌─────┐
r  │      │   │      │        │      │    │     │
──►│  Wi  ├──►│  K1  ├───+───►│  W1  ├───►│  G  ├─┐y
   │      │   │      │   │    │      │    │     │ │
   └──────┘   └──────┘   │    └──────┘    └─────┘ │
                         │                        │
                         │    ┌──────┐            │
                         │    │      │            │
                         └────┤  K2  ◄────────────┘
                              │      │
                              └──────┘
```
Where the returned controller `K` takes the measurement vector `[r; y]` (positive feedback), 
i.e., it includes all blocks `Wi, K1, K2, W1`.
If `match_dc = true`, `Wi` is automatically computed to make sure the static gain matches `Tref` exactly, otherwise `Wi` is set to `I`.
The `info` named tuple contains the feedforward filter for inspection (`info.K1 = K1*Wi`).


# Arguments:
- `G`: Plant model
- `Tref`: Reference model
- `γ`: Relative γ
- `ρ`: Design parameter, typically 1 < ρ < 3. Increase to emphasize model matching at the expense of robustness.
- `W1`: Pre-compensator for loop shaping.
- `Wo`: Output selction matrix. If there are more measurements than controlled variables, this matrix let's you select which measurements are to be controlled. 
- `kwargs`: Are sent to [`hinfsynthesize`](@ref).

Ref: Sec. 9.4.3 of Skogestad, "Multivariable Feedback Control: Analysis and Design".
The reference contains valuable pointers regarding gain-scheduling implementation of the designed controller as an observer with feedback from estimated states.
In order to get anti-windup protection when `W1` contains an integrator,
transform `W1` to self-conditioned Hanus form (using [`hanus`](@ref)) and implement the controller like this
```julia
W1h = hanus(W1)             # Perform outside loop

# Each iteration
us = filter(Ks, [r; y])     # filter inputs through info.Ks (filter is a fictive function that applies the transfer function)
u  = filter(W1h, [us; ua])  # filter us and u-actual (after input saturation) through W1h
ua = clamp(u, lower, upper) # Calculate ua for next iteration as the saturated value of u
```


# Example:
```julia
P = tf([1, 5], [1, 2, 10]) # Plant
W1 = tf(1,[1, 0]) |> ss    # Loop shaping controller

Tref = tf(1, [1, 1]) |> ss # Reference model

K1dof, γ1, info1 = glover_mcfarlane(ss(P), 1.1; W1)
K2dof, γ2, info2 = glover_mcfarlane_2dof(ss(P), Tref, 1.1, 1.1; W1)

G1 = feedback(P*K1dof)
G2 = info2.Gcl

bodeplot(info2.K1, w, lab="Feedforward filter")
plot([step(G1, 15), step(G2, 15), step(Tref, 15)], lab=["1-DOF" "2-DOF" "Tref"])
```
"""
function glover_mcfarlane_2dof(G::AbstractStateSpace{Continuous}, Tref::AbstractStateSpace{Continuous}, γ = 1.1, ρ = 1.1; W1=1, Wo = I, match_dc = true, kwargs...)
    γ > 1 || throw(ArgumentError("γ must be greater than 1"))
    ρ > 1 || throw(ArgumentError("ρ must be greater than 1"))
    Gs = G*W1
    As,Bs,Cs,Ds = ssdata(Gs)
    Ar,Br,Cr,Dr = ssdata(Tref)
    nr,nr = size(Ar)
    lr,mr = size(Dr)
    ns,ns = size(As)
    ls,ms = size(Ds)
    Rs = I + Ds*Ds'
    sRs = sqrt(Rs) # Example in Skogestad uses matlab sqrt which is elementwise
    Ss = lu!(I + Ds'Ds)
    A1 = (As - Bs*(Ss\Ds'Cs))
    R1 = Cs'*(Rs\Cs)
    Q1 = Bs*(Ss\Bs')
    Zs, _ = MatrixEquations.arec(A1', R1, Q1)

    A = cat(As, Ar, dims=(1,2))
    B1 = [
        zeros(ns,mr) ((Bs*Ds')+(Zs*Cs'))/sRs
        Br zeros(nr,ls)
    ]
    B2 = [Bs; zeros(nr,ms)]
    C1 = [zeros(ms,ns+nr); Cs zeros(ls,nr); ρ*Wo*Cs -ρ^2*Cr]
    C2 = [zeros(mr,ns+nr); Cs zeros(ls,nr)]
    D11 = [zeros(ms,mr+ls); zeros(ls,mr) sRs;-ρ^2*Dr ρ*Wo*sRs]
    D12 = [I(ms);Ds;ρ*Ds]
    D21 = [ρ*I(mr) zeros(mr,ls);zeros(ls,mr) sRs]
    D22 = [zeros(mr,ms); Ds]
    P = ss(A, B1, B2, C1, C2, D11, D12, D21, D22)
    Ks, γopt = hinfsynthesize(P, γrel = γ; kwargs...)
    
    u1 = 1:ms
    K1 = Ks[:, u1]
    if match_dc
        K2 = Ks[:, ms+1:end]
        sens = output_sensitivity(Gs, -K2)*Gs*K1 # eq. 9.89
        Wi = (Wo*dcgain(sens, 1e-6))\dcgain(Tref)
        K1 = K1*Wi
        Ks.B[:,u1] .= K1.B
        Ks.D[:,u1] .= K2.D
    else
        Wi = I
    end
    K = W1*Ks
    Gcl = feedback(G*K, ss(1), W1=1:ms, U1=ms+1:2ms, pos_feedback=true)
    info = (; Gcl, Ks, Gs, Wi, K1)
    K, γopt, info
end



"Spectral radius"
function ρ(X)
    e = eigvals(X)
    abs(e[end])
end


"""
    Wh = hanus(W)

Return `Wh` on Hanus form. `Wh` has twice the number of inputs, where the second half of the inputs are "actual inputs", e.g., potentially saturated. This is used to endow `W` with anti-windup protection.
`W` must have an invertable `D` matrix and be minimum phase.

Ref: Sec 9.4.5 of Skogestad, "Multivariable Feedback Control: Analysis and Design"
"""
function hanus(W)
    A,B,C,D = ssdata(W)
    nu = W.nu
    BD = B/D
    A2 = A - BD*C
    B2 = [0*I(nu) BD]
    D2 = [D 0*I(nu)]
    ss(A2, B2, C, D2)
end

"""
    extended_gangoffour(P, C)

Returns a single statespace system that maps 
- `w1` reference or measurement noise
- `w2` load disturbance
to
- `z1` control error
- `z2` control input
```
      z1          z2
      ▲  ┌─────┐  ▲      ┌─────┐
      │  │     │  │      │     │
w1──+─┴─►│  C  ├──┴───+─►│  P  ├─┐
    │    │     │      │  │     │ │
    │    └─────┘      │  └─────┘ │
    │                 w2         │
    └────────────────────────────┘
```

The returned system has the transfer-function matrix
```math
\\begin{bmatrix}
I \\\\ C
\\end{bmatrix} (I + PC)^{-1} \\begin{bmatrix}
I & P
\\end{bmatrix}
```

The gang of four can be plotted like so
```julia
Gcl = extended_gangoffour(G, C) # Form closed-loop system
bodeplot(Gcl, lab=["S" "CS" "PS" "T"], plotphase=false) |> display # Plot gang of four
```
Note, the last output of Gcl is the negative of the `CS` and `PS` transfer functions from `gangoffour2`.
See [`glover_mcfarlane`](@ref) for an extended example. See also [`ncfmargin`](@ref).
"""
function extended_gangoffour(P, C)
    ny,nu = size(P)
    S = feedback(ss(I(ny+nu), P.timeevol), [0*I(ny) P; -C 0*I(nu)], pos_feedback=true)
    Gcl = S + cat(0*I(ny), -I(nu), dims=(1,2))
    Gcl
end

"""
    m, ω = ncfmargin(P, K)

Normalized coprime factor margin, defined has the inverse H∞ norm of
```math
\\begin{bmatrix}
I \\\\ K
\\end{bmatrix} (I + PK)^{-1} \\begin{bmatrix}
I & P
\\end{bmatrix}
```
A margin ≥ 0.25-0.3 is a reasonable for robustness. 

If controller `K` stabilizes `P` with margin `m`, then `K` will also stabilize `P̃` if `nugap(P, P̃) < m`.

See also [`extended_gangoffour`](@ref), [`diskmargin`](@ref).
"""
function ncfmargin(P, K)
    Gcl = extended_gangoffour(P, K)
    im, w = hinfnorm2(Gcl)
    inv(im), w
end