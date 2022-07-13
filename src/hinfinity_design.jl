""" The implementation is primarily based on a paper by Glover and Doyle, and
the code is best read with this docment at hand. All references to equations
in docstrings and comments are to the version of the paper given below:

  @article{glover1988state,
    title={State-space formulae for all stabilizing controllers that satisfy an
           H-infinity norm bound and relations to relations to risk sensitivity},
    author={Glover, Keith and Doyle, John C},
    journal={Systems & control letters},
    volume={11},
    number={3},
    pages={167--172},
    year={1988},
    publisher={Citeseer}
  }
"""


"""
    flag = hinfassumptions(P::ExtendedStateSpace; verbose=true)

Check the assumptions for using the γ-iteration synthesis in Theorem 1.
"""
function hinfassumptions(P::ExtendedStateSpace; verbose = true)

    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(P)

    # Check assumption A1
    if !_stabilizable(A, B2)
        verbose && @warn(
            "The system A is not stabilizable through B2, ",
            "violation of assumption A1."
        )
        return false
    end
    if !_detectable(A, C2)
        verbose && @warn(
            "The system A is not detectable through C2, ",
            "violation of assumption A1."
        )
        return false
    end

    # Check assumption A2
    if rank(D12) < size(D12, 2)
        verbose &&
            @warn("The matrix D12 does not have full rank, ", "violation of assumption A2. The full control signal must have nonzero penalty at infinite frequency.")
        return false
    end
    if rank(D21) < size(D21, 1)
        verbose &&
            @warn("The matrix D21 does not have full rank, ", "violation of assumption A2. The whole measurement vector y must be corrupted by noise at infinite frequency.")
        return false
    end

    # Check assumption A5
    if rank(A - B2 * pinv(D12) * C1) < size(A, 1)
        verbose && @warn(
            "The matrix (A - B2*D12^-*C1) does not have full",
            "rank, violation of assumption A5."
        )
        return false
    end
    # Check assumption A6
    if rank(A - B1 * pinv(D21) * C2) < size(A, 1)
        verbose && @warn(
            "The matrix (A - B1*D21Pinv*C2) does not ",
            "have full rank, violation of assumption A6."
        )
        return false
    end

    # All assumptions have passed, and we may proceed with the synthesis
    verbose && println("All assumtions are satisfied!")
    return true
end

"""
    flag = _stabilizable(A::AbstractMatrix, B::AbstractMatrix)

Applies the Hautus lemma to check if the pair is stabilizable
"""
function _stabilizable(A::AbstractMatrix, B::AbstractMatrix)
    eigValsA = eigvals(A)
    for ii in eachindex(eigValsA)
        if real(eigValsA[ii]) >= 0
            if rank([eigValsA[ii] * I - A B]) != size(A, 1)
                return false
            end
        end
    end
    return true
end

"""
    flag = _detectable(A::AbstractMatrix, C::AbstractMatrix)

Applies the Hautus lemma to check if the pair is detectable
"""
function _detectable(A::AbstractMatrix, C::AbstractMatrix)
    eigValsA = eigvals(A)
    for ii = 1:length(eigValsA)
        if real(eigValsA[ii]) >= 0
            if rank([eigValsA[ii] * I - A; C]) != size(A, 1)
                return false
            end
        end
    end
    return true
end

"""
    K, γ, mats = hinfsynthesize(P::ExtendedStateSpace; gtol = 1e-4, interval = (0, 20), verbose = false, tolerance = 1.0e-10, γrel = 1.01, transform = true, ftype = Float64, check = true)

Computes an H-infinity optimal controller `K` for an extended plant `P` such that
||F_l(P, K)||∞ < γ (`lft(P, K)`) for the smallest possible γ given `P`. The routine is
known as the γ-iteration, and is based on the paper "State-space formulae for
all stabilizing controllers that satisfy an H∞-norm bound and relations to
risk sensitivity" by Glover and Doyle.


# Arguments:
- `gtol`: Tolerance for γ.
- `interval`: The starting interval for the bisection.
- `verbose`: Print progress?
- `tolerance`: For detecting eigenvalues on the imaginary axis.
- `γrel`: If `γrel > 1`, the optimal γ will be found by γ iteration after which a controller will be designed for `γ = γopt * γrel`. It is often a good idea to design a slightly suboptimal controller, both for numerical reasons, but also since the optimal controller may contain very fast dynamics. If `γrel → ∞`, the computed controller will approach the 𝑯₂ optimal controller. Getting a mix between 𝑯∞ and 𝑯₂ properties is another reason to choose `γrel > 1`.
- `transform`: Apply coordiante transform in order to tolerate a wider range or problem specifications.
- `ftype`: construct problem matrices in higher precision for increased numerical robustness. If the calculated controller achieves 
- `check`: Perform a post-design check of the γ value achieved by the calculated controller. A warning is issued if the achieved γ differs from the γ calculated during design. If this warning is issued, consider using a higher-precision number type like `ftype = BigFloat`.

See the example folder for example usage.
"""
function hinfsynthesize(
    P::ExtendedStateSpace{Continuous, T};
    gtol = 1e-4,
    interval = (0.0, 20.0),
    verbose = false,
    tolerance = sqrt(eps(T)),
    γrel = 1.01,
    transform = true,
    ftype = Float64,
    check = true, 
) where T
    Thigh = promote_type(T, ftype)
    hp = Thigh != T
    if hp
        bb(x) = Thigh.(x)
        mats = bb.(ssdata_e(P))
        Pa = ss(mats..., P.timeevol)
    else
        Pa = P
    end

    # Transform the system into a suitable form
    if transform
        P̄, Ltrans12, Rtrans12, Ltrans21, Rtrans21 = _transformp2pbar(Pa)
    else
        P̄, Ltrans12, Rtrans12, Ltrans21, Rtrans21 = Pa, I, I, I, I, I
    end

    # Run the γ iterations
    X∞Feasible, Y∞Feasible, F∞Feasible, H∞Feasible, γFeasible =
        _γiterations(P̄, interval, verbose, gtol, tolerance)


    if γFeasible !== nothing
        # Synthesize the controller and transform it back into the original coordinates

        if γrel > 1
            γFeasible *= γrel
            X∞Feasible, Y∞Feasible, F∞Feasible, H∞Feasible = _solvematrixequations(P̄, γFeasible)
        end

        K = _synthesizecontroller(
            P̄,
            X∞Feasible,
            Y∞Feasible,
            F∞Feasible,
            H∞Feasible,
            γFeasible,
            Ltrans12,
            Rtrans12,
            Ltrans21,
            Rtrans21,
        )

        # Return the controller, the optimal gain γ
        γ = γFeasible
    else
        # Return and empty controller, empty gain γ
        K = ss(0.0)
        γ = Inf
    end
    if hp
        bf(x) = T.(x)
        mats = bf.(ssdata(K))
        K = ss(mats..., K.timeevol)
    end
    if check
        γactual = hinfnorm2(lft(P, K))[1]::T
        diff = γ - γactual
        abs(diff) > 10gtol && @warn "Numerical problems encountered, returned γ is adjusted to the γ achieved by the computed controller (γ - γactual = $diff). Try solving the problem in higher precision by calling hinfsynthesize(...; ftype=BigFloat)"
        γFeasible = γactual
    end
    return K, γFeasible, (X=X∞Feasible, Y=Y∞Feasible, F=F∞Feasible, H=H∞Feasible, P̄)
end

"""
    K = _synthesizecontroller(P::ExtendedStateSpace, Xinf, Yinf, F, H, γ, Ltrans12, Rtrans12, Ltrans21, Rtrans21)

Syntheize a controller by operating on the scaled state-space description of the
system (i.e., the state-space realization of `P̄`) using the solutions from the
γ-iterations. The controller is synthesized in the coordinates of `P̄`, and then
transformed back to the original coordinates by the linear transformations
Ltrans12, Rtrans12, Ltrans21 and Rtrans21.
"""
function _synthesizecontroller(
    P::ExtendedStateSpace,
    Xinf,
    Yinf,
    F,
    H,
    γ::Number,
    Ltrans12,
    Rtrans12,
    Ltrans21,
    Rtrans21,
)

    A   = P.A
    B1  = P.B1
    B2  = P.B2
    C1  = P.C1
    C2  = P.C2
    D11 = P.D11
    D12 = P.D12
    D21 = P.D21
    D22 = P.D22

    γ² = γ * γ

    B = [B1 B2]
    C = [C1; C2]

    # Dimensionality
    P1 = size(C1, 1)
    P2 = size(C2, 1)
    M1 = size(B1, 2)
    M2 = size(B2, 2)

    # Equation (11)
    # F11 = F[1:(M1-P2), :]
    F12 = F[(M1-P2+1):M1, :]
    F2 = F[(M1+1):(M1+M2), :]

    # Equation (12)
    # H11 = H[:, 1:(P1-M2)]
    H12 = H[:, (P1-M2+1):P1]
    H2 = H[:, (P1+1):(P1+P2)]

    # Definition of D in the assumptions section
    D1111 = D11[1:(P1-M2), 1:(M1-P2)]
    D1112 = D11[1:(P1-M2), (M1-P2+1):M1]
    D1121 = D11[(P1-M2+1):P1, 1:(M1-P2)]
    D1122 = D11[(P1-M2+1):P1, (M1-P2+1):M1]

    # Equation 19
    J1 = (γ² * I - D1111 * D1111')
    D11hat = ((-D1121 * D1111') / J1) * D1112 - D1122

    # Equation 20
    D12hatD12hat = I - (D1121 / (γ² * I - D1111' * D1111)) * D1121'
    _assertrealandpsd(D12hatD12hat; msg = " in equation (20)")
    D12hat = cholesky(Hermitian(D12hatD12hat)).L

    # Equation 21
    D21hatD21hat = I - (D1112' / J1) * D1112
    _assertrealandpsd(D21hatD21hat; msg = " in equation (21)")
    D21hat = cholesky(Hermitian(D21hatD21hat)).U

    # Equation 27
    Zinv = (I - Yinf * Xinf / γ²)

    # Equation 22
    B2hat = (B2 + H12) * D12hat

    # Equation 23 (using the inverse of 27)
    C2hat = -D21hat * (C2 + F12) / Zinv

    # Equation 24
    B1hat = -H2 + (B2hat / D12hat) * D11hat

    # Equation 25 (using the inverse of 27)
    C1hat = F2 / Zinv + (D11hat / D21hat) * C2hat

    # Equation 26
    Ahat = A + H * C + (B2hat / D12hat) * C1hat

    Ac = Ahat

    B1c = B1hat * Ltrans21
    B2c = B2hat

    C1c = Rtrans12 * C1hat
    C2c = C2hat

    Bc = [B1c B2c]
    Cc = [C1c; C2c]

    D11c = Rtrans12 * D11hat * Ltrans21
    # TODO implement loop shift for any system not satisfying A4
    # D12c = D12hat * Ltrans21
    # D21c = Rtrans12 * D21hat
    # D22c = zeros(size(D11hat))
    # Dc = [D11c D12c; D21c D22c]
    # return ss(Ac, Bc[:, 1:P2], Cc[1:M2, :], Dc[1:M2, 1:P2])

    return ss(Ac, Bc[:, 1:P2], Cc[1:M2, :], D11c)
end

"""
    rqr(D, γ=1)
    
"Regularized" qr factorization. This struct represents \$(D'D + γI)\$ without forming \$D'D\$
Note: this does not support negative γ or \$(γI - D'D)\$
Supported operations: `\\,/,*`, i.e., it behaves also like a matrix (unlike the standard `QR` factorization object).
"""
struct rqr{T, PT, DGT}
    P::PT
    DG::DGT
    function rqr(D, γ=1)
        P = (D'D) + γ*I
        DG = qr([D; √(γ)I])
        new{eltype(P), typeof(P), typeof(DG)}(P, DG)
    end
end

Base.:\(d::rqr, b) = (d.DG.R\(adjoint(d.DG.R)\b))
Base.:/(b, d::rqr) = ((b/d.DG.R)/adjoint(d.DG.R))
Base.:*(d::rqr, b) = (d.P*b)
Base.:*(b, d::rqr) = (b*d.P)


"""
    _assertrealandpsd(A::AbstractMatrix, msg::AbstractString)

Check that a matrix is real and PSD - throw an error otherwise.
"""
function _assertrealandpsd(A::AbstractMatrix; msg = "")
    any(real(eigvals(A)) .<= 0) && error(string("The matrix", msg, " is not PSD."))
    any(imag(eigvals(A)) .!= 0) && error(string("The matrix", msg, " is not real."))
end

"""
    flag =  _checkfeasibility(Xinf, Yinf, γ, tolerance, iteration; verbose=true)

Check the feasibility of the computed solutions Xinf, Yinf and the algebraic
Riccatti equations, return true if the solution is valid, and false otherwise.
"""
function _checkfeasibility(
    Xinf::AbstractMatrix,
    Yinf::AbstractMatrix,
    γ::Number,
    tolerance::Number,
    iteration::Number;
    verbose = true,
)

    minXev = minimum(real.(eigvals(Xinf)))
    minYev = minimum(real.(eigvals(Yinf)))
    specrad = maximum(abs.(eigvals(Xinf * Yinf))) / (γ^2)

    if verbose
        iteration == 1 && println("iteration, γ")
        println(iteration, " ", γ)
    end

    # Failed test, eigenvalues of Xinf must be positive real
    minXev < -tolerance && return false
    # Failed test, eigenvalues of Yinf must be positive real
    minYev < -tolerance && return false
    # Failed test, spectral radius of XY must be greater than γ squared
    specrad > 1 && return false
    return true
end

"""
    solution = _solvehamiltonianare(H)

Solves a hamiltonian Alebraic Riccati equation using the Schur-decomposition,
for additional details, see

  @article{laub1979schur,
    title={A Schur method for solving algebraic Riccati equations},
    author={Laub, Alan},
    journal={IEEE Transactions on automatic control},
    volume={24},
    number={6},
    pages={913--921},
    year={1979},
    publisher={IEEE}
  }
"""
function _solvehamiltonianare(H::AbstractMatrix{T}) where T
    S = schur(H)
    So = ordschur(S, real.(S.values) .< 0)::Schur{T, typeof(H)}
    Z::Matrix{T} = So.Z
    (m, n) = size(Z)
    U11 = Z[1:div(m, 2), 1:div(n, 2)]
    U21 = Z[div(m, 2)+1:m, 1:div(n, 2)]

    return U21 / (U11), So.values # Note: if pinv is used, BigFloats may fail
end

"""
    solution = _solvematrixequations(P::ExtendedStateSpace, γ::Number)

Solves the dual matrix equations in the γ-iterations (equations 7-12 in Doyle).
"""
function _solvematrixequations(P::ExtendedStateSpace, γ::Number)
    A = P.A
    B1 = P.B1
    B2 = P.B2
    C1 = P.C1
    C2 = P.C2
    D11 = P.D11
    D12 = P.D12
    D21 = P.D21
    D22 = P.D22

    T = float(eltype(A))

    P1 = size(C1, 1)
    P2 = size(C2, 1)
    M1 = size(B1, 2)
    M2 = size(B2, 2)

    γ² = γ * γ

    B = [B1 B2]
    C = [C1; C2]

    # Equation (7)
    D1dot = [D11 D12]
    R = [Matrix(-γ²*I(M1)) zeros(T, M1, M2); zeros(T, M2, M1) zeros(T, M2, M2)] + D1dot' * D1dot |> svd # Explicit call to Matrix constructor required for type stability https://github.com/JuliaLang/julia/issues/44408

    # Equation (8)
    Ddot1 = [D11; D21]
    Rbar::Matrix{T} = [Matrix(-γ²*I(P1)) zeros(T, P1, P2); zeros(T, P2, P1) zeros(T, P2, P2)] + Ddot1 * Ddot1' 

    # Form hamiltonian for X and Y, equation (9) and (10)
    HX::Matrix{T} = [A zeros(T, size(A)); -C1'*C1 -A'] - ([B; -C1' * D1dot]) * (R \ [D1dot' * C1 B'])
    HY::Matrix{T} = [A' zeros(T, size(A)); -B1*B1' -A] - ([C'; -B1 * Ddot1']) * (svd(Rbar) \ [Ddot1 * B1' C]) # cat with Adjoint is type unstable in julia v1.7 but not in v1.9, hence the typeassert

    # Solve matrix equations
    Xinf, vx = _solvehamiltonianare(HX)
    Yinf, vy = _solvehamiltonianare(HY)

    # Equation (11)
    F = -(R \ (D1dot' * C1 + B' * Xinf))

    # Equation (12)
    H = -(B1 * Ddot1' + Yinf * C') / Rbar

    return Xinf, Yinf, F, H, vx, vy
end

"""
    X∞, Y∞, F, H, γ = _γiterations(P, interval, verbose, gtol, tolerance)

Rune the complete set of γ-iterations over a specified search interval with a
set number of iterations. It is possible to break the algorithm if the number
of iterations become too large. This should perhaps be tken care of by
specifying an interval and tolerance for γ. In addition, the algorithm simply
terminates without a solution if the maximum possible γ on the defined
interval is infeasible.
"""
function _γiterations(
    P::ExtendedStateSpace,
    interval::Tuple,
    verbose::Bool,
    gtol, 
    tolerance::Number,
)

    T = typeof(P.A)
    ET = eltype(T)
    XinfFeasible, YinfFeasible, FinfFeasible, HinfFeasible, gammFeasible =
        T(undef,0,0), T(undef,0,0), T(undef,0,0), T(undef,0,0), nothing

    gl, gu = ET.(interval)
    gl = max(ET(1e-3), gl)
    iters = ceil(Int, log2((gu-gl+1e-16)/gtol))  

    for iteration = 1:iters
        γ = sqrt(gu*gl)
        # Solve the matrix equations
        Xinf, Yinf, F, H, vx, vy = _solvematrixequations(P, γ)

        # Check Feasibility

        isFeasible = all(abs.(real.(vx)) .> tolerance) && all(abs.(real.(vy)) .> tolerance) &&
            _checkfeasibility(Xinf, Yinf, γ, tolerance, iteration; verbose = verbose)

        if isFeasible
            XinfFeasible = Xinf
            YinfFeasible = Yinf
            FinfFeasible = F
            HinfFeasible = H
            gammFeasible = γ
            gu = γ
        else
            gl = γ
        end
    end
    return XinfFeasible, YinfFeasible, FinfFeasible, HinfFeasible, gammFeasible
end

"""
    Pbar, Ltrans12, Rtrans12, Ltrans21, Rtrans21 = _transformP2Pbar(P::ExtendedStateSpace)

Transform the original system P to a new system Pbar, in which D12bar = [0; I]
and D21bar = [0 I] in order to satisfy the feasibility assumption A3 (see Doyle)
"""
function _transformp2pbar(P::ExtendedStateSpace)

    # Compute the transformation
    Ltrans12, Rtrans12 = _scalematrix(P.D12)
    Ltrans21, Rtrans21 = _scalematrix(P.D21)

    # Transform the system
    Abar = P.A
    B1bar = P.B1 * Rtrans21
    B2bar = P.B2 * Rtrans12
    C1bar = Ltrans12 * P.C1
    C2bar = Ltrans21 * P.C2
    D11bar = Ltrans12 * P.D11 * Rtrans21
    D12bar = Ltrans12 * P.D12 * Rtrans12
    D21bar = Ltrans21 * P.D21 * Rtrans21
    D22bar = Ltrans21 * P.D22 * Rtrans12
    Pbar = ss(Abar, B1bar, B2bar, C1bar, C2bar, D11bar, D12bar, D21bar, D22bar)

    return Pbar, Ltrans12, Rtrans12, Ltrans21, Rtrans21
end

"""
    Tl, Tr = _scalematrix(A::AbstractMatrix; method::String)

Find a left and right transform of A such that Tl*A*Tr = [I, 0], or
Tl*A*Tr = [I; 0], depending on the dimensionality of A.
"""
function _scalematrix(A::AbstractMatrix; method = :QR)
    # Check the rank condition
    if (minimum(size(A)) > 0)
        if rank(A) != minimum(size(A))
            error("Cannot scale the system, assumption A2 is violated. This typically means that the control input is not penalized by D12 or that there is missing measurement noise in D21.")
        end
    else
        error("Cannot scale the system, minimum size of A must begreater than 0")
    end

    # Perform scaling with the cosen method
    if method === :QR
        return _coordinatetransformqr(A)
    elseif method === :SVD
        return _coordinatetransformsvd(A)
    else
        error("The method $method is not supported, use 'QR' or 'SVD' instad.")
    end
end

"""
    Tl, Tr =  _computeCoordinateTransformQR(A::AbstractMatrix)

Use the QR decomposition to find a transformaiton [Tl, Tr] such that
Tl*A*Tr becomes [0;I], [0 I] or I depending on the dimensionality of A.
"""
function _coordinatetransformqr(A::AbstractMatrix)
    m, n = size(A)
    if m == n
        # Square matrix with full rank
        Q, R = qr(A)
        LeftTransform = Q'
        RightTransform = inv(R)
    elseif m > n
        # Rectangular matrix with rank = n
        Q, R = qr(A)
        LeftTransform = [Q[:, (n+1):(m)]'; Q[:, 1:n]']
        RightTransform = inv(R[1:n, 1:n])
    elseif m < n
        # Rectangular matrix with rank = m
        Q, R = qr(A')
        LeftTransform = inv(R[1:m, 1:m])'
        RightTransform = [Q[:, m+1:n] Q[:, 1:m]]
    end
    return LeftTransform, RightTransform
end

"""
    Tl, Tr =  _computeCoordinateTransformSVD(A::AbstractMatrix)

Use the SVD to find a transformaiton [Tl, Tr] such that
Tl*A*Tr becomes [0;I], [0 I] or I depending on the dimensionality of A.
"""
function _coordinatetransformsvd(A::AbstractMatrix)
    m, n = size(A)
    if m == n
        # Square matrix with full rank
        U, S, V = svd(A)
        LeftTransform = inv(Diagonal(S)) * U'
        RightTransform = V
    elseif m > n
        # Rectangular matrix with rank = n
        U, S, V = svd(A; full = true)
        LeftTransform = [U[:, n+1:m]'; inv(Diagonal(S)) * U[:, 1:n]']
        RightTransform = V
    elseif m < n
        # Rectangular matrix with rank = m
        U, S, V = svd(A; full = true)
        LeftTransform = inv(Diagonal(S)) * U'
        RightTransform = [V[:, m+1:n] V[:, 1:m]]
    end
    return LeftTransform, RightTransform
end

"""
    P = hinfpartition(G, WS, WU, WT)

Transform a SISO or MIMO system G, with weighting functions WS, WU, WT into
an LFT with an isolated controller, and write the resulting system, P(s),
on a state-space form. Valid inputs for G are transfer functions (with dynamics,
can be both MIMO and SISO, both in tf and ss forms). Valid inputs for the
weighting functions are empty arrays, numbers (static gains), and `LTISystem`s.

Note, `system_mapping(P)` is equal to `-G`.
"""
function hinfpartition(G, WS, WU, WT)
    WS isa LTISystem && common_timeevol(G,WS)
    WU isa LTISystem && common_timeevol(G,WU)
    WT isa LTISystem && common_timeevol(G,WT)
    te = G.timeevol
    # # Convert the systems into state-space form
    Ag, Bg, Cg, Dg = _input2ss(G)
    Asw, Bsw, Csw, Dsw = _input2ss(WS)
    Auw, Buw, Cuw, Duw = _input2ss(WU)
    Atw, Btw, Ctw, Dtw = _input2ss(WT)

    G = ss(Ag, Bg, Cg, Dg)
    WS = ss(Asw, Bsw, Csw, Dsw)
    WU = ss(Auw, Buw, Cuw, Duw)
    WT = ss(Atw, Btw, Ctw, Dtw)

    # Check that the system is realizable
    if size(Cg, 1) != size(Btw, 2) && size(Btw, 2) != 0
        if ControlSystems.issiso(WT)
            WT = ControlSystems.append(fill(WT, G.ny)...)
            return hinfpartition(G, WS, WU, WT)
        end
        println([size(Cg, 1), size(Btw, 2)])
        throw(DimensionMismatch(
            "You must have the same number of inputs to WT as there are outputs",
        ))
    end
    if size(Cg, 1) != size(Bsw, 2) && size(Bsw, 2) != 0
        if ControlSystems.issiso(WS)
            WS = ControlSystems.append(fill(WS, G.ny)...)
            return hinfpartition(G, WS, WU, WT)
        end
        println([size(Cg, 1), size(Bsw, 2)])
        throw(DimensionMismatch(
            "You must have the same number of inputs to WS as there are outputs",
        ))
    end
    if size(Bg, 2) != size(Buw, 2) && size(Buw, 2) != 0
        if ControlSystems.issiso(WU)
            WU = ControlSystems.append(fill(WU, G.nu)...)
            return hinfpartition(G, WS, WU, WT)
        end
        println([size(Bg, 2), size(Buw, 2)])
        throw(DimensionMismatch(
            "You must have the same number of inputs to WU as there are controls",
        ))
    end
    if (
        size(Ag, 1) == 0 ||
        size(Ag, 2) == 0 ||
        size(Bg, 1) == 0 ||
        size(Bg, 2) == 0 ||
        size(Cg, 1) == 0 ||
        size(Cg, 2) == 0 ||
        size(Dg, 1) == 0 ||
        size(Dg, 2) == 0
    )
        throw(DimensionMismatch(
            "Expansion of systems dimensionless A,B,C or D is not yet supported",
        ))
    end

    # Form A
    (mAg, nAg) = size(Ag)
    (mAsw, nAsw) = size(Asw)
    (mAuw, nAuw) = size(Auw)
    (mAtw, nAtw) = size(Atw)

    if size(Bsw, 1) != 0 && size(Bsw, 2) != 0
        BswCg = Bsw * Cg
    else
        BswCg = zeros(mAsw, nAg)
    end
    if size(Btw, 1) != 0 && size(Btw, 2) != 0
        BtwCg = Btw * Cg
    else
        BtwCg = zeros(mAtw, nAg)
    end

    @debug ([(mAg, nAg), (mAsw, nAsw), (mAuw, nAuw), (mAtw, nAtw)])
    A = [
        Ag zeros(mAg, nAsw) zeros(mAg, nAuw) zeros(mAg, nAtw)
        -BswCg Asw zeros(mAsw, nAuw) zeros(mAsw, nAtw)
        zeros(mAuw, nAg) zeros(mAuw, nAsw) Auw zeros(mAuw, nAtw)
        BtwCg zeros(mAtw, nAsw) zeros(mAtw, nAuw) Atw
    ]

    if size(Buw, 2) == 0
        Buw = zeros(0, size(Bg, 2))
    end

    (mBg, nBg) = size(Bg)
    (mBsw, nBsw) = size(Bsw)
    (mBuw, nBuw) = size(Buw)
    (mBtw, nBtw) = size(Btw)

    Bw = [zeros(mBg, nBsw); Bsw; zeros(mBuw, nBsw); zeros(mAtw, nBsw)]
    Bu = [Bg; zeros(mBsw, nBg); Buw; zeros(mAtw, nBg)]

    (mCg, nCg) = size(Cg)
    (mCsw, nCsw) = size(Csw)
    (mCuw, nCuw) = size(Cuw)
    (mCtw, nCtw) = size(Ctw)

    if size(Dsw, 1) != 0 && size(Dsw, 1) != 0
        DswCg = Dsw * Cg
    else
        DswCg = zeros(0, nAg)
    end
    if size(Dtw, 1) != 0 && size(Dtw, 1) != 0
        DtwCg = Dtw * Cg
    else
        DtwCg = zeros(0, nAg)
    end

    Cz = [
        -DswCg Csw zeros(mCsw, nAuw) zeros(mCsw, nAtw)
        zeros(mCuw, nAg) zeros(mCuw, nAsw) Cuw zeros(mCuw, nAtw)
        DtwCg zeros(mCtw, nAsw) zeros(mCtw, nAuw) Ctw
    ]
    Cy = [-Cg zeros(mCg, nAsw) zeros(mCg, nAuw) zeros(mCg, nAtw)]


    if size(Duw, 2) == 0
        Duw = zeros(0, size(Bg, 2))
    end

    (mDg, nDg) = size(Dg)
    (mDsw, nDsw) = size(Dsw)
    (mDuw, nDuw) = size(Duw)
    (mDtw, nDtw) = size(Dtw)

    Dzw = [Dsw; zeros(mDuw, nDsw); zeros(mDtw, nDsw)]
    Dzu = [zeros(mDsw, nDuw); Duw; zeros(mDtw, nDuw)]
    Dyw = Matrix{Float64}(I, mCg, nDuw)
    Dyu = -Dg

    P = ss(A, Bw, Bu, Cz, Cy, Dzw, Dzu, Dyw, Dyu, te)

end

# function hinfpartition2(G, WS, WU, WT)
#     # # Convert the systems into state-space form
#     Ag, Bg, Cg, Dg = _input2ss(G)
#     Asw, Bsw, Csw, Dsw = _input2ss(WS)
#     Auw, Buw, Cuw, Duw = _input2ss(WU)
#     Atw, Btw, Ctw, Dtw = _input2ss(WT)

#     G = ss(Ag, Bg, Cg, Dg)
#     WS = ss(Asw, Bsw, Csw, Dsw)
#     WU = ss(Auw, Buw, Cuw, Duw)
#     WT = ss(Atw, Btw, Ctw, Dtw)

#     # Check that the system is realizable
#     if size(Cg, 1) != size(Btw, 2) && size(Btw, 2) != 0
#         if ControlSystems.issiso(WT)
#             WT = ControlSystems.append(fill(WT, G.ny)...)
#             return hinfpartition(G, WS, WU, WT)
#         end
#         println([size(Cg, 1), size(Btw, 2)])
#         throw(DimensionMismatch(
#             "You must have the same number of inputs to WT as there are outputs",
#         ))
#     end
#     if size(Cg, 1) != size(Bsw, 2) && size(Bsw, 2) != 0
#         if ControlSystems.issiso(WS)
#             WS = ControlSystems.append(fill(WS, G.ny)...)
#             return hinfpartition(G, WS, WU, WT)
#         end
#         println([size(Cg, 1), size(Bsw, 2)])
#         throw(DimensionMismatch(
#             "You must have the same number of inputs to WS as there are outputs",
#         ))
#     end
#     if size(Bg, 2) != size(Buw, 2) && size(Buw, 2) != 0
#         if ControlSystems.issiso(WU)
#             WU = ControlSystems.append(fill(WU, G.nu)...)
#             return hinfpartition(G, WS, WU, WT)
#         end
#         println([size(Bg, 2), size(Buw, 2)])
#         throw(DimensionMismatch(
#             "You must have the same number of inputs to WU as there are controls",
#         ))
#     end
#     if (
#         size(Ag, 1) == 0 ||
#         size(Ag, 2) == 0 ||
#         size(Bg, 1) == 0 ||
#         size(Bg, 2) == 0 ||
#         size(Cg, 1) == 0 ||
#         size(Cg, 2) == 0 ||
#         size(Dg, 1) == 0 ||
#         size(Dg, 2) == 0
#     )
#         throw(DimensionMismatch(
#             "Expansion of systems dimensionless A,B,C or D is not yet supported",
#         ))
#     end


#     ny,nu = G.ny, G.nu
#     Iny = I(ny)
#     perm1 = [Iny; zeros(ny+nu,ny); Iny]
#     perm2 = [zeros(ny,nu) ; I(nu) ; zeros(2*ny,nu)] + [-Iny ; zeros(nu,ny) ; Iny ; -Iny] * G
#     W = ControlSystems.blockdiag(WS, WU, WT, ss(Iny)) 
#     P = W * [perm1 perm2] 
#     nz = P.ny - ny
#     nw = P.nu - nu
#     zi = 1:nz
#     yi = nz+1:P.ny
#     wi = 1:nw
#     ui = nw+1:P.nu
#     ss(P.A, P.B[:, wi], P.B[:, ui], P.C[zi, :], P.C[yi, :], 
#         P.D[zi, wi], P.D[zi, ui], P.D[yi, wi], P.D[yi, ui], G.timeevol)
# end

function ControlSystems.blockdiag(systems::AbstractStateSpace...)
    ST = promote_type(typeof.(systems)...)
    timeevol = common_timeevol(systems...)
    A = ControlSystems.blockdiag(s.A for s in systems)
    B = ControlSystems.blockdiag(_remove_empty_cols(s.B) for s in systems)
    C = ControlSystems.blockdiag(s.C for s in systems)
    D = ControlSystems.blockdiag(_remove_empty_cols(s.D) for s in systems)
    return ST(A, B, C, D, timeevol)
end

_remove_empty_cols(x) = size(x,2) == 0 ? zeros(size(x,1), 1) : x

ControlSystems.numeric_type(n::Number) = typeof(n)
ControlSystems.numeric_type(::Any) = Float64
ControlSystems.numeric_type(::AbstractVector{T}) where T = numeric_type(T)

"""
    convert_input_to_ss(H)

Helper function used for type conversion in hinfpartition()
"""
function _input2ss(H)
    T = ControlSystems.numeric_type(H) |> float
    if isa(H, LTISystem)
        if isa(H, TransferFunction)
            Hss = ss(H)
        else
            Hss = H
        end
        Ah, Bh, Ch, Dh = Hss.A, Hss.B, Hss.C, Hss.D
    elseif isa(H, Number)
        Ah, Bh, Ch, Dh = zeros(T, 0, 0), zeros(T, 0, 1), zeros(T, 1, 0), H * ones(T, 1, 1)
    else # e.g. H = nothing
        Ah, Bh, Ch, Dh = zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(0, 0)
    end
    return Ah, Bh, Ch, Dh
end

"""
    hinfsignals(P::ExtendedStateSpace, G::LTISystem, C::LTISystem)

Use the extended state-space model, a plant and the found controller to extract
the closed loop transfer functions.

- `Pcl : w → z` : From input to the weighted functions
- `S   : w → e` : From input to error
- `CS  : w → u` : From input to control
- `T   : w → y` : From input to output
"""
function hinfsignals(P::ExtendedStateSpace, G::LTISystem, C::LTISystem)

    common_timeevol(P,G,C)
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(P)
    _, _, Cg, Dg = ssdata(ss(G))
    Ac, Bc, Cc, Dc = ssdata(ss(C))

    # Precompute the inverse
    M = inv(I - Dc * D22)
    A11 = A + B2 * M * Dc * C2
    A12 = B2 * M * Cc
    A21 = Bc * C2 + Bc * D22 * M * Dc * C2
    A22 = Ac + Bc * D22 * M * Cc
    B11 = B1 + B2 * M * Dc * D21
    B21 = Bc * D21 + Bc * D22 * M * Dc * D21

    C_w2z_11 = C1 + D12 * M * Dc * C2
    C_w2z_12 = D12 * M * Cc
    D_w2z = D11 + D12 * M * Dc * D21

    Pw2z = ss([A11 A12; A21 A22], [B11; B21], [C_w2z_11 C_w2z_12], D_w2z, timeevol(P))

    C_w2e_11 = C2 + D22 * M * Dc * C2
    C_w2e_12 = D22 * M * Cc
    D_w2e = D21 + D22 * M * Dc * D21

    Pw2e = ss([A11 A12; A21 A22], [B11; B21], [C_w2e_11 C_w2e_12], D_w2e, timeevol(P))

    C_w2u_11 = M * Dc * C2
    C_w2u_12 = M * Cc
    D_w2u = M * Dc * D21

    Pw2u = ss([A11 A12; A21 A22], [B11; B21], [C_w2u_11 C_w2u_12], D_w2u, timeevol(P))

    Abar = [A11 A12; A21 A22]
    Bbar = [B11; B21]
    C_w2y_11 = M * Dc * C2
    C_w2y_12 = M * Cc
    D_w2y = M * Dc * D21
    Cbar1 = [Cg zeros(size(Cg, 1), size(Abar, 2) - size(Cg, 2))]
    Cbar2 = [Dg * C_w2y_11 Dg * C_w2y_12]
    Dbar = Dg * M * Dc * D21

    Pw2y = ss(Abar, Bbar, Cbar1 + Cbar2, Dbar, timeevol(P))

    Pcl = Pw2z # Verified to be equal to julia> lft(P, C) == Pcl -> true
    S = Pw2e
    CS = Pw2u
    T = Pw2y

    return Pcl, S, CS, T
end


"""
    bilineard2c(Ad::AbstractArray, Bd::AbstractArray, Cd::AbstractArray, Dd::AbstractArray, Ts::Number; tolerance=1e-12)

Balanced Bilinear transformation in State-Space. This method computes a
continuous time equivalent of a discrete time system, such that

    G_c(z) = z2s[G_d(z)]

in a manner which accomplishes the following
  (i)   Preserves the infinity L-infinity norm over the transformation
  (ii)  Finds a system which balances B and C, in the sense that ||B||_2=||C||_2
  (iii) Satisfies G_d(z) = s2z[z2s[G_d(z)]] for some map s2z[]
"""
function bilineard2c(
    Ad::AbstractArray,
    Bd::AbstractArray,
    Cd::AbstractArray,
    Dd::AbstractArray,
    Ts::Number;
    tolerance = 1e-12,
)

    Id = Matrix{Float64}(I, size(Ad, 1), size(Ad, 2))

    Pd = Ad - Id
    Qd = Ad + Id
    ialpha = 2 / Ts #Should be this, but the nyquist frequency doesnt add up unless
    ialpha = 1 / Ts

    Ac = ialpha * (Pd / Qd)
    Bc = Qd \ Bd
    Cc = 2 * ialpha * (Cd / Qd)
    Dc = Dd - (Cd / Qd) * Bd

    # Scaling for improved numerical stability
    σB = maximum(svd(Bc).S)
    σC = maximum(svd(Cc).S)
    if σB > tolerance && σC > tolerance
        λd = sqrt(σB / σC)
    else
        λd = 1
        error(
            "Warning, the problem is poorly cnditioned. Consider an alternate discretization scheme.",
        )
    end
    Bc /= λd
    Cc *= λd

    return Ac, Bc, Cc, Dc
end

"""
    bilineard2c(sys::StateSpace)

Applies a Balanced Bilinear transformation to continuous-time statespace object
"""
function bilineard2c(sys::StateSpace)
    Ad, Bd, Cd, Dd = ssdata(sys)
    Ts = sys.Ts

    if Ts <= 0
        error("Error, the input must be a discrete time system.")
    end

    Ac, Bc, Cc, Dc = bilineard2c(Ad, Bd, Cd, Dd, Ts)
    return ss(Ac, Bc, Cc, Dc)
end

"""
    bilineard2c(sys::ExtendedStateSpace)

Applies a Balanced Bilinear transformation to continuous-time extended statespace object
"""
function bilineard2c(sys::ExtendedStateSpace)
    Ad = sys.A
    Bd = sys.B
    Cd = sys.C
    Dd = sys.D
    Ts = sys.Ts

    m1 = size(sys.B1, 2)
    m2 = size(sys.B2, 2)
    p1 = size(sys.C1, 1)
    p2 = size(sys.C2, 1)

    if Ts <= 0
        error("Error, the input must be a discrete time system.")
    end

    Ac, Bc, Cc, Dc = bilineard2c(Ad, Bd, Cd, Dd, Ts)

    A = Ac
    B1 = Bc[:, 1:m1]
    B2 = Bc[:, (m1+1):(m1+m2)]
    C1 = Cc[1:p1, :]
    C2 = Cc[(p1+1):(p1+p2), :]
    D11 = Dc[1:p1, 1:m1]
    D12 = Dc[1:p1, (m1+1):(m1+m2)]
    D21 = Dc[(p1+1):(p1+p2), 1:m1]
    D22 = Dc[(p1+1):(p1+p2), (m1+1):(m1+m2)]

    return ss(A, B1, B2, C1, C2, D11, D12, D21, D22)
end

"""
    bilinearc2d(Ac::AbstractArray, Bc::AbstractArray, Cc::AbstractArray, Dc::AbstractArray, Ts::Number; tolerance=1e-12)

Balanced Bilinear transformation in State-Space. This method computes a
discrete time equivalent of a continuous-time system, such that

\$G_d(z) = s2z[G_c(s)]\$

in a manner which accomplishes the following
  (i)   Preserves the infinity L-infinity norm over the transformation
  (ii)  Finds a system which balances B and C, in the sense that \$||B||_2=||C||_2\$
  (iii) Satisfies \$G_c(s) = z2s[s2z[G_c(s)]]\$ for some map z2s[]
"""
function bilinearc2d(
    Ac::AbstractArray,
    Bc::AbstractArray,
    Cc::AbstractArray,
    Dc::AbstractArray,
    Ts::Number;
    tolerance = 1e-12,
)

    Id = Matrix{Float64}(I, size(Ac, 1), size(Ac, 2))
    alpha = Ts / 2 #Should be this, but the nyquist frequency doesnt add up
    alpha = Ts

    # Check that the bilinear tranformation is possible
    if minimum(svd(Id - alpha * Ac).S) < 1e-12
        error(
            "The transformation is extremely poorly conditioned, with min(svd(Id - alpha * Ac).S) < 1e-12. Consider an alternate discretization scheme.",
        )
    end

    PP = Id - alpha * Ac
    QQ = Id + alpha * Ac

    Ad = PP \ QQ
    Bd = (PP \ Bc)
    Cd = 2 * alpha * (Cc / PP)

    # Scaling for improved numerical stability
    σB = maximum(svd(Bd).S)
    σC = maximum(svd(Cd).S)
    if σB > tolerance && σC > tolerance
        λc = sqrt(σB / σC)
    else
        λc = 1
        error(
            "Warning, the problem is poorly cnditioned. Consider an alternate discretization scheme.",
        )
    end

    Bd /= λc
    Cd *= λc

    Dd = alpha * Cc / PP * Bc + Dc
    return Ad, Bd, Cd, Dd, Ts
end

"""
    bilinearc2d(sys::StateSpace, Ts::Number)

Applies a Balanced Bilinear transformation to a discrete-time statespace object
"""
function bilinearc2d(sys::StateSpace{Continuous}, Ts::Number)
    Ac, Bc, Cc, Dc = ssdata(sys)

    if Ts <= 0
        throw(ArgumentError("Error, the the discretization time Ts must be positive."))
    end

    Ad, Bd, Cd, Dd = bilinearc2d(Ac, Bc, Cc, Dc, Ts)
    return ss(Ad, Bd, Cd, Dd, Ts)
end

"""
    bilinearc2d(sys::ExtendedStateSpace, Ts::Number)

Applies a Balanced Bilinear transformation to a discrete-time extended statespace object
"""
function bilinearc2d(sys::ExtendedStateSpace{Continuous}, Ts::Number)
    Ac = sys.A
    Bc = sys.B
    Cc = sys.C
    Dc = sys.D

    m1 = size(sys.B1, 2)
    m2 = size(sys.B2, 2)
    p1 = size(sys.C1, 1)
    p2 = size(sys.C2, 1)

    if Ts <= 0
        error("Error, the the discretization time Ts must be positive.")
    end

    Ad, Bd, Cd, Dd = bilinearc2d(Ac, Bc, Cc, Dc, Ts)

    A = Ad
    B1 = Bd[:, 1:m1]
    B2 = Bd[:, (m1+1):(m1+m2)]
    C1 = Cd[1:p1, :]
    C2 = Cd[(p1+1):(p1+p2), :]
    D11 = Dd[1:p1, 1:m1]
    D12 = Dd[1:p1, (m1+1):(m1+m2)]
    D21 = Dd[(p1+1):(p1+p2), 1:m1]
    D22 = Dd[(p1+1):(p1+p2), (m1+1):(m1+m2)]

    return ss(A, B1, B2, C1, C2, D11, D12, D21, D22, Ts)
end

function fudge_inv(s::AbstractStateSpace, ε = 1e-3)
    s = deepcopy(s)
    s.D .+= ε*I(size(s.D,1))
    inv(s)
end


function hinfgrad(sys; rtolinf=1e-8, kwargs...)
    hn, ω = hinfnorm2(sys; rtolinf, kwargs...)
    ∇ = hinfgrad(sys, hn, ω)
    (∇..., hn, ω)
end

"""
    ∇A, ∇B, ∇C, ∇D, hn, ω = hinfgrad(sys; rtolinf=1e-8, kwargs...)
    ∇A, ∇B, ∇C, ∇D        = hinfgrad(sys, hn, ω)

Compute the gradient of the H∞ norm w.r.t. the statespace matrices `A,B,C,D`.
If only a system is provided, the norm `hn` and the peak frequency `ω` are automatically calculated. `kwargs` are sent to [`hinfnorm2`](@ref).
Note, the default tolerance to which the norm is calculated is set smaller than default for [`hinfnorm2`](@ref), gradients will be discontinuous with any non-finite tolerance, and sensitive optimization algorithms may require even tighter tolerance.

In cases where the maximum singular value is reached at more than one frequency, a random frequency is used.

If the system is unstable, the gradients are `NaN`. Strategies to find an initial stabilizing controllers are outlined in Apkarian and D. Noll, "Nonsmooth H∞ Synthesis" in IEEE Transactions on Automatic Control.

An `rrule` for ChainRules is defined using this function, so `hn` is differentiable with any AD package that derives its rules from ChainRules (only applies to the `hn` return value, not `ω`).
"""
function hinfgrad(sys, hn, ω)
    A,B,C,D = ssdata(sys)
    isfinite(hn) ||
        return fill(NaN, size(A)), fill(NaN, size(B)), fill(NaN, size(C)), fill(NaN, size(D))
        
    if !isfinite(ω)
        Γ⁻¹B = zeros(size(B))
        CΓ⁻¹ = zeros(size(C))
    else
        # Γ = complex(0, ω)*I - A
        Γ = ω*I - A
        Γf = lu(Γ)
        Γ⁻¹B = Γf\B
        CΓ⁻¹ = C/Γf
    end
    U,_,V = svd(C*Γ⁻¹B + D)
    UV = U[:,1]*V[:,1]'
    ∇A = CΓ⁻¹'*UV*Γ⁻¹B'
    ∇B = CΓ⁻¹'*UV
    ∇C = UV*Γ⁻¹B'
    ∇D = UV
    ∇A, ∇B, ∇C, ∇D
end

function ChainRulesCore.rrule(hinfnorm::H, sys::T; rtolinf=1e-8, kwargs...) where {T <: AbstractStateSpace,
    H <: Union{typeof(hinfnorm), typeof(hinfnorm2)}
}
    hn, ω = hinfnorm(sys; rtolinf, kwargs...)
    function hinfnorm2_pullback((Δhn, Δw))
        # @show Δhn, Δw
        Δw isa ZeroTangent || error("adjoint for the frequency returned by $H not implemented")
        ∇A, ∇B, ∇C, ∇D = hinfgrad(sys, hn, ω)
        ∇D .*= Δhn
        ∇C .*= √(Δhn)
        ∇B .*= √(Δhn)
        senssys = Tangent{T}(; A=∇A, B=∇B, C=∇C, D=∇D, timeevol=NoTangent())
        NoTangent(), senssys
    end
    return (hn, ω), hinfnorm2_pullback
end
