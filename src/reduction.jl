using ControlSystems: ssdata
# TODO:  Implement Reducing unstable linear control systems via real Schur transformation.
# By temporarily removing unstable dynamics, reduce, add dynamics back.

# TODO: consider implementing methods on slides 158-161 https://cscproxy.mpi-magdeburg.mpg.de/mpcsc/benner/talks/lecture-MOR.pdf
# this requires a generalization of the method below to take a function from 
# sys -> P, Q

"""
    frequency_weighted_reduction(G, Wo, Wi)

Find Gr such that ||Wₒ(G-Gr)Wᵢ||∞ is minimized.
For a realtive reduction, set Wo = inv(G) and Wi = I.

Ref: Robust and Optimal Control ch. 7.2
"""
function frequency_weighted_reduction(G, Wo, Wi, r)
    A,B,C,D = ssdata(G)

    sys = Wo*G*Wi
    n = nstates(G)
    
    if Wi == 1 || Wi == I
        P = gram(G, :c)
    else
        P0 = gram(sys, :c)
        m = size(P0, 1) - n
        P = [I(n) zeros(n, m)] * P0 * [I(n); zeros(m, n)]
    end
    if Wo == 1 || Wo == I
        Q = gram(G, :o)
    else
        Q0 = gram(sys, :o)
        m = size(Q0, 1) - n
        Q = [I(n) zeros(n, m)] * Q0 * [I(n); zeros(m, n)]
    end

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
    R = factorize(I - D22*Dk)
    R̃ = factorize(I - Dk*D22)
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
    controller_reduction(P, K, r, out=false)

Minimize    ||(K-Kᵣ) W||∞ if out=false
            ||W (K-Kᵣ)||∞ if out=true
See Robust and Optimal Control Ch 19.1
out indicates if the weight will be applied as output or input weight.
"""
function controller_reduction(P, K, r, out=false)
    W = controller_reduction_weight(P, K)
    if out
        frequency_weighted_reduction(K, W, 1, r)
    else
        frequency_weighted_reduction(K, 1, W, r)
    end
end

"""
    hsvd(sys::AbstractStateSpace{Continuous})

Return the Hankel singular values of `sys`, computed as the eigenvalues of `QP`
Where `Q` and `P` are the Gramians of `sys`.
"""
function hsvd(sys::AbstractStateSpace{Continuous})
    P = gram(sys, :c)
    Q = gram(sys, :o)
    e = eigvals(Q * P)
    sqrt.(e)
end

# slide 189 https://cscproxy.mpi-magdeburg.mpg.de/mpcsc/benner/talks/lecture-MOR.pdf
# This implementation works, but results in a complex-valued system.
# function model_reduction_irka(sys::AbstractStateSpace{Continuous}, r; tol = 1e-6)
#     A,B,C,D = ssdata(sys)
#     all(iszero, D) || error("Nonzero D not supported.")
#     sysr0, _ = baltrunc(sys, n=r)
#     Ah, Bh, Ch = real.(ssdata(sysr0))
#     e, T = eigen(Ah)
#     eold = 100e
#     k = 0
#     while maximum(abs.(e .- eold) ./ abs.(e)) > tol && k < 1000
#         Bt = T'Bh
#         Ct = Ch*T
#         Vv = map(1:r) do i
#             (-e[i]*I - A)\(B*Bt[i, :])
#         end
#         V = reduce(hcat, Vv)
#         Wv = map(1:r) do i
#             (-e[i]*I - A')\(C'*Ct[:,i])
#         end
#         W = reduce(hcat, Wv)
#         V = orthonormal(V) # not sure about those
#         W = orthonormal(W)
#         Ah = (W'V)\W'A*V
#         Bh = (W'V)\W'B
#         Ch = C*V
#         eold = e
#         k += 1
#         e, T = eigen(Ah)
#     end
#     sysr = ss(Ah, Bh, Ch, 0)

# end

# function orthonormal(X)
#     # U,S,V = svd(X, full=true)
#     # R = U*V'
#     Matrix(qr(X).Q) # |> real # this works sometimes
# end