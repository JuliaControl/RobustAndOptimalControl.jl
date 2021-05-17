"""
    δ(N=32)

Create an uncertain element of `N` uniformly distributed samples ∈ [-1, 1]
"""
δ(N=32) = StaticParticles(N, Uniform(-1, 1))

function ControlSystems.lft(M11::AbstractMatrix, M12::AbstractMatrix, M21::AbstractMatrix, M22::AbstractMatrix, d::AbstractMatrix)
    M11 + M12 * d*((I - M22*d)\M21)
end

wass(a, b) = wasserstein(a,b,2)
l2(a,b) = sum(abs2.(a.particles .- b.particles))

wass(a::AbstractArray, b::AbstractArray) = sum(wass(a,b) for (a,b) in zip(a,b)) # shouldn't be allowed to permute all particles individually?
l2(a::AbstractArray, b::AbstractArray) = sum(l2(a,b) for (a,b) in zip(a,b)) # shouldn't be allowed to permute all particles individually?


function pars2lft(M11, p, δ)
    @unpack M12, M21, M22 = p
    Phat = lft(M11, M12, M21, M22, δ)
end

struct LFT
    M11
    M12
    M21
    M22
    d
end

Base.Matrix(l::LFT) = [l.M11 l.M12; l.M21 l.M22]

function Base.show(io::IO, l::LFT)
    display(Matrix(l))
end

function ControlSystems.lft(l::LFT)
    @unpack M11, M12, M21, M22, d = l
    lft(M11, M12, M21, M22, d)
end

function ControlSystems.ss(l::LFT, P0::AbstractStateSpace)
    P = lft(l)
    A = P[:, 1:P0.nx]
    B = P[:, P0.nx+1:end]
    ss(A, B, P0.C, P0.D)
end

find_lft(sys::StateSpace{<:Any, <:StaticParticles{<:Any, N}}, n_uncertain::Int)  where N = find_lft(sys, [δ(N) for _ in 1:n_uncertain], wass)

"""
    l = find_lft(sys::StateSpace{<:Any, <:StaticParticles{<:Any, N}}, δ) where N

Given an systems `sys` with uncertain coefficients in the form of `StaticParticles`, find a lower linear fractional transformation `M` such that `lft(M, δ) ≈ sys`. 

`δ` can be either the source of uncertainty in `sys`, i.e., a vector of the unique uncertain parameters that were used to create `sys`. These should be constructed as uniform randomly distributed particles for most robust-control theory to be applicable. 
`δ` can also be an integer, in which case a numer of `δ` sources of uncertainty are automatically created. This could be used for order reduction if the number of uncertainty sources in `sys` is large.

Note, uncertainty in `sys` is only supported in `A` and `B`, `C` and `D` must be deterministic.

Returns `l::LFT` that internaly contains all four blocks of `M` as well as `δ`. Call `ss(l,sys)` do obtain `lft(M, δ) ≈ sys`.

Call `Matrix(l)` to obtain `M = [M11 M12; M21 M22]`
"""
function find_lft(sys::StateSpace{<:Any, <:StaticParticles{<:Any, N}}, delta, dist::F = l2) where {N, F}
    n_uncertain = length(delta)
    P = [sys.A sys.B]
    M11 = mean.(P)
    δ = Diagonal(delta)
    function loss(p)
        Phat = pars2lft(M11, p, δ)
        dist(P, Phat) #+ 1e-7*sum(abs, p)
    end
    p0 = ComponentVector(
        M12 = 0.00001randn(sys.nx, n_uncertain),
        M21 = 0.00001randn(n_uncertain, sys.nx+sys.nu),
        M22 = 0.00001randn(n_uncertain, n_uncertain),
    )
    res = Optim.optimize(
        loss,
        p0,
        LBFGS(),
        Optim.Options(
            store_trace       = true,
            show_trace        = true,
            show_every        = 5,
            iterations        = 40000,
            allow_f_increases = false,
            time_limit        = 45,
            x_tol             = 1e-8,
            f_reltol          = 0,
            g_tol             = 1e-8,
            f_calls_limit     = 0,
            g_calls_limit     = 0,
        ),
        # autodiff = :forward
    )
    p = res.minimizer
    LFT(M11, p.M12, p.M21, p.M22, δ), res

end


