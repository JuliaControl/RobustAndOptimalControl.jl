
"""
    K, Cl = h2synthesize(P::ExtendedStateSpace, γ = nothing)

Synthesize H₂-optimal controller K and calculate the closed-loop transfer function from `w` to `z`.
Ref: Cha. 14.5 in Robust and Optimal Control.

If `γ = nothing`, use the formulas for H₂ in Ch 14.5. If γ is a large value, the H∞ formulas are used. As γ → ∞, these two are equivalent. The h∞ formulas do a coordinate transfromation that handles slightly more general systems so if you run into an error, it might be worth trying setting γ to something large, e.g., 1000.
"""
function h2synthesize(P::ExtendedStateSpace, γ = nothing)

    if γ === nothing
        X2, Y2, F2, L2 = _solvematrixequations2(P)
        Â2 = P.A + P.B2*F2 + L2*P.C2
        K = ss(Â2, -L2, F2, 0)
        return K, lft(ss(P), K)
    end
    
    P̄, Ltrans12, Rtrans12, Ltrans21, Rtrans21 = _transformp2pbar(P)
    X2, Y2, F2, L2 = _solvematrixequations(P̄, γ)

    _checkfeasibility(
        X2,
        Y2,
        γ,
        1e-3,
        10;
        verbose = false,
    ) || throw(DomainError("Not feasible for γ = $γ"))


    K = _synthesizecontroller(
        P̄,
        X2,
        Y2,
        F2,
        L2,
        γ,
        Ltrans12,
        Rtrans12,
        Ltrans21,
        Rtrans21,
    )
    Cl = lft(ss(P), K)
    K, Cl
end

function _solvematrixequations2(P::ExtendedStateSpace)

    A = P.A
    B1 = P.B1
    B2 = P.B2
    C1 = P.C1
    C2 = P.C2
    D11 = P.D11
    D12 = P.D12
    D21 = P.D21
    D22 = P.D22

    P1 = size(C1, 1)
    P2 = size(C2, 1)
    M1 = size(B1, 2)
    M2 = size(B2, 2)

    HX = [A zeros(size(A)); -C1'*C1 -A'] - [B2; -C1'*D12] * [D12'*C1 B2']
    HY = [A' zeros(size(A)); -B1*B1' -A] - [C2'; -B1*D21'] * [D21*B1' C2]

    # Solve matrix equations
    X2 = _solvehamiltonianare(HX)
    Y2 = _solvehamiltonianare(HY)

    F2 = -(B2'X2 + D12'C1)

    L2 = -(B1 * D21' + Y2 * C2')

    return X2, Y2, F2, L2
end