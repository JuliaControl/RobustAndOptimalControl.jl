function blockdiagonalize(A::AbstractMatrix)
    E = eigen(A, sortby=eigsortby)
    Db,Vb = cdf2rdf(E)
    Db,Vb,E
end

eigsortby(λ::Real) = λ
eigsortby(λ::Complex) = (abs(imag(λ)),real(λ))

function complex_indices(A::Matrix) # assumes A on block diagonal form
    findall(diag(A, -1) .!= 0)
end

function real_indices(A::Matrix) # assumes A on block diagonal form
    size(A,1) == 1 && return [1]
    setdiff(findall(diag(A, -1) .== 0), complex_indices(A).+1)
end

function complex_indices(D::AbstractVector)
    complex_eigs = imag.(D) .!= 0
    findall(complex_eigs)
end

function cdf2rdf(E::Eigen)
    # Implementation inspired by scipy https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cdf2rdf.html
    # with the licence https://github.com/scipy/scipy/blob/v1.6.3/LICENSE.txt
    D,V = E
    n = length(D)

    # get indices for each first pair of complex eigenvalues
    complex_inds = complex_indices(D)

    # eigvals are sorted so conjugate pairs are next to each other
    j = complex_inds[1:2:end]
    k = complex_inds[2:2:end]

    # put real parts on diagonal
    Db = zeros(n, n)
    Db[diagind(Db)] .= real(D)

    # compute eigenvectors for real block diagonal eigenvalues
    U = zeros(eltype(D), n, n)
    U[diagind(U)] .= 1.0

    # transform complex eigvals to real blockdiag form
    for (k,j) in zip(k,j)
        Db[j, k] = imag(D[j]) # put imaginary parts in blocks 
        Db[k, j] = imag(D[k])

        U[j, j] = 0.5im
        U[j, k] = 0.5
        U[k, j] = -0.5im
        U[k, k] = 0.5
    end
    Vb = real(V*U)

    return Db, Vb
end

"""
    sysm, T, E = modal_form(sys; C1 = false)

Bring `sys` to modal form.

The modal form is characterized by being tridiagonal with the real values of eigenvalues of `A` on the main diagonal and the complex parts on the first sub and super diagonals. `T` is the similarity transform applied to the system such that 
```julia
sysm ≈ similarity_transform(sys, T)
```

If `C1`, then an additional convention for SISO systems is used, that the `C`-matrix coefficient of real eigenvalues is 1. If `C1 = false`, the `B` and `C` coefficients are chosen in a balanced fashion.

`E` is an eigen factorization of `A`.

See also [`hess_form`](@ref) and [`schur_form`](@ref)
"""
function modal_form(sys; C1 = false)
    Ab,T,E = blockdiagonalize(sys.A)
    # Calling similarity_transform looks like a detour, but this implementation allows modal_form to work with any AbstractStateSpace which implements a custom method for similarity transform
    sysm = similarity_transform(sys, T)
    sysm.A .= Ab # sysm.A should already be Ab after similarity_transform, but Ab has less numerical noise
    if ControlSystems.issiso(sysm)
        # This enforces a convention: the C matrix entry for the first component in each mode is positive. This allows SISO systems on modal form to be interpolated in a meaningful way by interpolating their coefficients. 
        # Ref: "New Metrics Between Rational Spectra and their Connection to Optimal Transport" , Bagge Carlson,  Chitre
        ci = complex_indices(sysm.A)
        flips = ones(sysm.nx)
        for i in ci
            if sysm.C[1, i] < 0
                flips[i] = -1
                flips[i .+ 1] = -1
            end
        end
        ri = real_indices(sysm.A)
        for i in ri
            c = sysm.C[1, i]
            if C1
                if c != 0
                    flips[i] /= c
                end
            else
                b = sysm.B[i, 1]
                flips[i] *= sqrt(abs(b))/(sqrt(abs(c)) + eps(b))
            end
        end
        T2 = diagm(flips)
        sysm = similarity_transform(sysm, T2)
        T = T*T2
        sysm.A .= Ab # Ab unchanged by diagonal T
    end
    sysm, T, E
end

"""
    sysm, T, SF = schur_form(sys)

Bring `sys` to Schur form.

The Schur form is characterized by `A` being Schur with the real values of eigenvalues of `A` on the main diagonal. `T` is the similarity transform applied to the system such that 
```julia
sysm ≈ similarity_transform(sys, T)
```
`SF` is the Schur-factorization of `A`.

See also [`modal_form`](@ref) and [`hess_form`](@ref)
"""
function schur_form(sys)
    SF = schur(sys.A)
    A = SF.T
    B = SF.Z'*sys.B
    C = sys.C*SF.Z
    ss(A,B,C,sys.D, sys.timeevol), SF.Z, SF
end

"""
    sysm, T, HF = hess_form(sys)

Bring `sys` to Hessenberg form form.

The Hessenberg form is characterized by `A` having upper Hessenberg structure. `T` is the similarity transform applied to the system such that 
```julia
sysm ≈ similarity_transform(sys, T)
```
`HF` is the Hessenberg-factorization of `A`.

See also [`modal_form`](@ref) and [`schur_form`](@ref)
"""
function hess_form(sys)
    F = hessenberg(sys.A)
    Q = Matrix(F.Q)
    A = F.H
    B = Q'sys.B 
    C = sys.C*Q
    D = sys.D
    ss(A,B,C,D, sys.timeevol), Q, F
end