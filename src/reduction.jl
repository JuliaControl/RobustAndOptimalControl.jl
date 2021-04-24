using ControlSystems: ssdata
"""
    frequency_weighted_reduction(G, Wo, Wi)

Find Gr such that ||Wₒ(G-Gr)Wᵢ||∞ is minimized.
For a realtive reduction, set Wo = inv(G) and Wi = I.
"""
function frequency_weighted_reduction(G, Wo, Wi, r)
    A,B,C,D = ssdata(G)
    @assert all(iszero, D) " not supported, but can be fixed see Ch. 7.2 Robust and Optimal Control"
    # Ao,Bo,Co,Do = ssdata(Wo)
    # Ai,Bi,Ci,Di = ssdata(Wi)
    # balreal(G)

    sys = Wo*G*Wi

    P0 = gram(sys, :c)
    Q0 = gram(sys, :o)

    n = nstates(G)
    m = size(P0, 1) - n

    P = [I(n) zeros(n, m)] * P0 * [I(n); zeros(m, n)]
    Q = [I(n) zeros(n, m)] * Q0 * [I(n); zeros(m, n)]


    L = cholesky(Hermitian((Q+Q')./2), check=false)
    issuccess(L) || @warn("Balanced realization failed: Observability grammian not positive definite, system needs to be observable. Result may be inaccurate.")

    Q1 = L.U
    U,Σ,V = svd(Q1*P*Q1')
    Σ .= sqrt.(Σ)
    Σ1 = diagm(0 => sqrt.(Σ))
    T = Σ1\(U'Q1)
    A,B,C,D = T*A/T, T*B, C/T, D
    ss(A[1:r, 1:r], B[1:r, :], C[:, 1:r], D, G.timeevol)
end


"
Lemma 19.1 See Robust and Optimal Control Ch 19.1
"
function controller_reduction_weight(P::ExtendedStateSpace, K)
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(P)
    Ak,Bk,Ck,Dk = ssdata(K)
    R = I - D22*Dk
    R̃ = I - Dk*D22
    Aw = [
            A+B2*Dk*(R\C2) B2*(R̃\Ck)
            Bk*(R\C2)      Ak+Bk*D22*(R̃\Ck)
    ]
    Bw = [
        B2/R̃
        Bk*D22/R̃
    ]
    Cw = [R\C2   R\D22*Ck]
    Dw = D22/R̃
    ss(Aw, Bw, Cw, Dw)
end

"""
    controller_reduction(P, K, r)

Minimize ||(K-Kᵣ) W||∞
See Robust and Optimal Control Ch 19.1
"""
function controller_reduction(P, K, r)
    W = controller_reduction_weight(P, K)
    frequency_weighted_reduction(K, ss(I(noutputs(K))), W, r)
end