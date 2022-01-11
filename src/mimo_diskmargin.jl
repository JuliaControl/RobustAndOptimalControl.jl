using ControlSystems, MonteCarloMeasurements
using Optim

# NOTE: for complex structured perturbations the typical upper bound algorithm appears to be quite tight. For real perturbations maybe the interval method has something to offer?

# function Base.:∈(a::Real, b::Complex{<:AbstractParticles})
#     mi, ma = pextrema(b.re)
#     mi <= a <= ma || return false
#     mi, ma = pextrema(b.im)
#     mi <= a <= ma
# end


# # any0det(D::Matrix{<:Complex{<:Interval}}) = 0 ∈ det(D)

# function any0det(D::Matrix{<:Complex{<:Interval}})
#     # CD = complex.(IntervalArithmetic.sup.(real.(D)), IntervalArithmetic.sup.(imag.(D)))
#     # dmid = det(CD)
#     # real(dmid) < 0 && return true
#     0 ∈ det(D)
# end

# """
#     bisect_a(P, K, w; W = (:), Z = (:), tol = 0.001, N = 32, upper = true, δ = δc, σ=-1)
#     bisect_a(M, D, w; W = (:), Z = (:), tol = 0.001, N = 32, upper = true, δ = δc, σ=-1)

# EXPERIMENTAL AND SUBJECT TO BUGS, BREAKAGE AND REMOVAL

# For each frequency in `w`, find the largest `a` such that the loop with uncertainty elements of norm no greater than `a`, located at the inputs `W` and outputs `Z` of `P`, is stable.

# By default, a conservative lower bound on the disk margin is returned. The conservatism comes from multiple factors
# - Each uncertainty element is represented as an IntervalBox in the complex plane rather than a ball with radius `a`.
# - IntervalArithmetic is inherently conservative due to the "dependency problem", i.e., `x-x ≠ 0`.

# If `upper = true`, a Monte-Carlo approach is used to find an upper bound of the disk margin, this will be optimistic. If you calculate both and they are tight, you have a good indication of the true disk margin.

# # Arguments:
# - `P`: Plant
# - `K`: Controller
# - `w`: Frequency vector
# - `W`: The inputs at which to insert a perturbation, defaults to all.
# - `Z`: The outputs at which to insert a perturbation, defaults to all.
# - `tol`: The tolerance in the bisection and the resolution of the resulting `a`.
# - `N`: The number of samples for the upper bound estimation. 
# - `upper`: Calculate upper or lower bound?
# - `δ = δc` for complex perturbations and `δ = δr` for real perturbations.
# """
# function bisect_a(M0::AbstractArray, D::AbstractMatrix; tol=1e-3, kwargs...)
#     @views map(axes(M0, 3)) do i
#         M = @views M0[:,:,i]
#         au = opnorm(M)
#         al = eltype(D) <: Complex ? ρ(M) : 0
#         au-al < tol && return al
#         iters = ceil(Int, log2((au-al)/tol))
#         local a
#         for j = 1:iters
#             a = (au+al)/2
#             if any0det(I - (1/a*M)*D) # 0 ∈ det(I - (a*M0[:,:,i])*D)
#                 al = a
#             else
#                 au = a
#             end
#         end
#         1/al
#     end
# end

# function bisect_a(args...; tol=1e-3, kwargs...)
#     M0, D = get_M(args...; kwargs...)
#     bisect_a(M0, D; tol=1e-3, kwargs...)
# end

# """
#     M,D = get_M(P, K, w; W = (:), Z = (:), N = 32, upper = false, δ = δc, σ=-1)

# Return the frequency response of `M` in the `M-Δ` formulation that arises when individual, complex/real perturbations are introduced on inputs `W` and outputs `Z` (defaults to all).

# # Arguments:
# - `P`: System
# - `K`: Controller
# - `w`: Frequency vector
# - `W`: Input indices that ar perturbed
# - `Z`: Output indices that ar perturbed
# - `N`: Number of samples for the upper bound computation
# - `upper`: Indicate whether an upper or lower bound is to be computed
# - `δ = δc` for complex perturbations and `δ = δr` for real perturbations.
# """
# function get_M(P::LTISystem, K::LTISystem, w::AbstractVector; W = (:), Z = (:), N = 32, upper=false, δ = δc, σ = -1)
#     Z1 = W2 = Z == (:) ? (1:P.ny) : Z
#     W1 = Z2 = W == (:) ? (1:P.nu) : W
#     ny,nu = length(Z1), length(W1)
#     D = Δ(ny+nu, δ)
#     if upper
#         D = rand(D, N)
#     else
#         D = Diagonal([Interval(d) for d in diag(D)])
#     end
#     if isempty(Z1) # No outputs
#         L = (K*P)[W1, W1]
#     elseif isempty(W1) # No inputs
#         L = (P*K)[Z1, Z1]
#     else
#         L = [ss(zeros(P.ny, P.ny), P.timeevol) P;-K ss(zeros(K.ny, K.ny), P.timeevol)]
#         L = L[[Z1; P.ny .+ Z2], [W1; P.nu .+ W2]]
#     end
#     n = L.ny
#     X = ss(kron([(1+σ)/2 -1;1 -1], I(n)), L.timeevol)
#     M = starprod(X,L)
#     # M = feedback(P, K; W2, Z2, Z1, W1) # TODO: this is probably not correct
#     M0 = freqresp(M, w)
#     M0 = permutedims(M0, (2,3,1))
#     M0, D
# end

# function get_M(L::LTISystem, w::AbstractVector; N = 32, upper=false, δ = δc, σ = -1)
#     ny,nu = size(L)
#     D = Δ(ny+nu, δ)
#     if upper
#         D = rand(D, N)
#     else
#         D = Diagonal([Interval(d) for d in diag(D)])
#     end
#     n = L.ny
#     X = ss(kron([(1+σ)/2 -1;1 -1], I(n)), L.timeevol)
#     M = starprod(X,L)
#     # M = feedback(P, K; W2, Z2, Z1, W1) # TODO: this is probably not correct
#     M0 = freqresp(M, w)
#     M0 = permutedims(M0, (2,3,1))
#     M0, D
# end

# function muloss(M, d)
#     D = Diagonal(d)
#     opnorm(D\M*D)
# end

#=
If there's a single block only, we have
μ = maximum(abs(e) for e in eigvals(M) if e is real) if the block is real*I
μ = ρ(M) if the block is complex*I
μ = opnorm(M) if the block is full complex


=#

"""
    μ = structured_singular_value(M; tol=1e-4, scalings=false, dynamic=false)

Compute (an upper bound of) the structured singular value μ for diagonal Δ of complex perturbations (other structures of Δ are not yet supported).
`M` is assumed to be an (n × n × N_freq) array or a matrix.

We currently don't have any methods to compute a lower bound, but if all perturbations are complex the spectral radius `ρ(M)` is always a lower bound (usually not a good one).

If `scalings = true`, return also a `n × nf` matrix `Dm` with the diagonal scalings `D` such that
```
D = Diagonal(Dm[:, i])
σ̄(D\\M[:,:,i]*D)
```
is minimized.

If `dynamic = true`, the perturbations are assumed to be time-varying `Δ(t)`. In this case, the same scaling is used for all frequencies and the returned `D` if `scalings=true` is a vector `d` such that `D = Diagonal(d)`.
"""
function structured_singular_value(M::AbstractArray{T}; tol=1e-4, scalings=false, dynamic=false) where T
    Ms1 = similar(M[:,:,1])
    Ms2 = similar(M[:,:,1])
    function muloss(M, d)
        D = Diagonal(d)
        mul!(Ms1, M, D)
        ldiv!(Ms2, D, Ms1)
        opnorm(Ms2)
    end
    n = size(M, 2)
    d0 = ones(real(T), n)
    opts = Optim.Options(
        store_trace       = false,
        show_trace        = false,
        show_every        = 10,
        iterations        = 1000,
        allow_f_increases = false,
        time_limit        = 100,
        x_tol             = 0,
        f_abstol          = 0,
        g_tol             = tol,
        f_calls_limit     = 0,
        g_calls_limit     = 0,
    )
    if dynamic
        res = Optim.optimize(
            d->maximum(i->muloss(@view(M[:,:,i]),d), axes(M, 3)),
            d0,
            # i == 1 ? ParticleSwarm() : BFGS(alphaguess = LineSearches.InitialStatic(alpha=0.9), linesearch = LineSearches.HagerZhang()),
            n == 1 ? BFGS() : ParticleSwarm(),
            opts,
        )
        d0 = res.minimizer
        return scalings ? (res.minimum, d0) : res.minimum
    else
        if scalings
            Dm = Matrix{Float64}(undef, n, size(M,3))
        end
        mu = map(axes(M, 3)) do i
            @views M0 = M[:,:,i]
            res = Optim.optimize(
                d->muloss(M0,d),
                d0,
                # i == 1 ? ParticleSwarm() : BFGS(alphaguess = LineSearches.InitialStatic(alpha=0.9), linesearch = LineSearches.HagerZhang()),
                n == 1 ? BFGS() : i == 1 ? ParticleSwarm() : NelderMead(), # Initialize using Particle Swarm
                opts,
            )
            d0 = res.minimizer # update initial guess
            if scalings
                Dm[:, i] .= d0
            end
            res.minimum
        end
        if scalings
            (M isa AbstractMatrix ? mu[] : mu), Dm
        else
            M isa AbstractMatrix ? mu[] : mu
        end
    end
end

"""
    robstab(M0::UncertainSS, w=exp10.(LinRange(-3, 3, 1500)); kwargs...)

Return the robust stability margin of an uncertain model, defined as the inverse of the structured singular value.
Currently, only diagonal complex perturbations supported.
"""
robstab(M0::UncertainSS, args...; kwargs...) = 1/norm(structured_singular_value(M0, args...; kwargs...), Inf)

"""
    structured_singular_value(M0::UncertainSS, [w::AbstractVector]; kwargs...)

- `w`: Frequency vector, if none is provided, the maximum μ over a brid 1e-3 : 1e3 will be returned.
"""
function structured_singular_value(M0::LTISystem, w::AbstractVector; kwargs...)
    if M0 isa UncertainSS
        if length(M0.Δ) == 1 && M0.Δ[] isa δDyn
            return sigma(M0.M, w)[1][:, 1]
        end
        all(d isa δ{<:Complex} for d in M0.Δ) || error("Structured singular value only supported for diagonal complex perturbations")
        M0 = M0.M # not it's own method due to method ambiguity misery
    end
    M = freqresp(M0, w)
    M = permutedims(M, (2,3,1))
    μ = structured_singular_value(M; kwargs...)
end

function structured_singular_value(M0::LTISystem; kwargs...)
    w = exp10.(LinRange(-3, 3, 1500))
    μ = structured_singular_value(M0, w; kwargs...)
    maximum(μ)
end

"""
    sim_diskmargin(L, σ::Real, w::AbstractVector)
    sim_diskmargin(L, σ::Real = 0)

Simultaneuous diskmargin at the outputs of `L`. 
Uses should consider using [`diskmargin`](@ref).
"""
function sim_diskmargin(L::LTISystem,σ::Real,w::AbstractVector)
    # S̄ = S+(σ-1)/2*I = lft([(1+σ)/2 -1;1 -1], L)
    n = L.ny
    X = ss(kron([(1+σ)/2 -1;1 -1], I(n)), L.timeevol)
    M = starprod(X,L)
    M0 = permutedims(freqresp(M, w), (2,3,1))
    mu = structured_singular_value(M0)
    imu = inv.(structured_singular_value(M0))
    simultaneous = [Diskmargin(imu, σ; ω0 = w, L) for (imu, w) in zip(imu,w)]
end

"""
    sim_diskmargin(P::LTISystem, C::LTISystem, σ::Real = 0)

Simultaneuous diskmargin at both outputs and inputs of `P`.
Ref: "An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet
https://arxiv.org/abs/2003.04771
"""
function sim_diskmargin(P::LTISystem, C::LTISystem, σ::Real=0)
    L = [ss(zeros(P.ny, P.ny)) P;-C ss(zeros(C.ny, C.ny))]
    sim_diskmargin(L,σ)
end

"""
    sim_diskmargin(L, σ::Real = 0)

Return the smallest simultaneous diskmargin over the grid 1e-3:1e3
"""
function sim_diskmargin(L, σ::Real=0)
    m = sim_diskmargin(L, σ, LinRange(-3, 3, 500))
    m = argmin(d->d.α, m)
end


"""
    diskmargin(P::LTISystem, C::LTISystem, σ, w::AbstractVector, args...; kwargs...)

Simultaneuous diskmargin at outputs, inputs and input/output simultaneously of `P`. 
Returns a named tuple with the fields `input, output, simultaneous_input, simultaneous_output, simultaneous` where `input` and `output` represent loop-at-a-time margins, `simultaneous_input` is the margin for simultaneous perturbations on all inputs and `simultaneous` is the margin for perturbations on all inputs and outputs simultaneously.
"""
function diskmargin(P::LTISystem, C::LTISystem, σ, w::AbstractVector, args...; kwargs...)
    L = C*P
    issiso(L) && return diskmargin(L, σ, w, args...)
    input = loop_diskmargin(L, σ, w)
    simultaneous_input = sim_diskmargin(L,σ,w)
    L = P*C
    output = loop_diskmargin(L, σ, w)
    simultaneous_output = sim_diskmargin(L,σ,w)
    te = P.timeevol
    L = [ss(zeros(P.ny, P.ny), te) P;-C ss(zeros(C.ny, C.ny), te)]
    simultaneous = sim_diskmargin(L,σ,w)
    (; input, output, simultaneous_input, simultaneous_output, simultaneous)
end


##

"""
    loop_diskmargin(L, args...; kwargs...)

Calculate the loop-at-a-time diskmargin for each output of `L`.

See also [`diskmargin`](@ref), [`sim_diskmargin`](@ref).
Ref: "An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet
https://arxiv.org/abs/2003.04771
"""
function loop_diskmargin(L::LTISystem, args...; kwargs...)
    dms = map(1:L.ny) do i
        open_L = broken_feedback(L, i)
        diskmargin(open_L, args...; kwargs...)
    end
end

"""
    loop_diskmargin(P, C, args...; kwargs...)

Calculate the loop-at-a-time diskmargin for each output and input of `P`.
See also [`diskmargin`](@ref), [`sim_diskmargin`](@ref).
Ref: "An Introduction to Disk Margins", Peter Seiler, Andrew Packard, and Pascal Gahinet
https://arxiv.org/abs/2003.04771
"""
function loop_diskmargin(P::LTISystem,C::LTISystem,args...; kwargs...)
    input = loop_diskmargin(C*P, args...; kwargs...)
    output = loop_diskmargin(P*C, args...; kwargs...)
    (; input, output)
end

"""
    broken_feedback(L, i)
Closes all loops in square MIMO system `L` except for loops `i`.
Forms L1 in fig 14. of "An Introduction to Disk Margins" https://arxiv.org/abs/2003.04771
"""
function broken_feedback(L::LTISystem, i)
    ny, nu = size(L)
    ny == nu || throw(ArgumentError("Only square loop-transfer functions supported"))
    connection_inds = setdiff(1:ny, i)
    i isa AbstractVector || (i = [i])
    open_L = feedback(
        L,
        ss(I(ny-1), L.timeevol), # open one loop
        U1 = connection_inds,
        Y1 = connection_inds,
        Z1 = i,
        W1 = i,
    )
    @assert issiso(open_L)
    open_L
end



"""
    loop_scaling(M0::Matrix, tol = 0.0001)

Find the optimal diagonal scaling matrix `D` such that `D\\M0*D` has a minimized condition number.
Applicable to square `M0` only. See also [`structured_singular_value`](@ref) with option `dynamic=true`.
Use [`loop_scale`](@ref) to find and apply the scaling to a loop-transfer function.
"""
function loop_scaling(M::Matrix, tol=1e-4)

    Ms1 = similar(M)
    Ms2 = similar(M)
    function scaleloss(d)
        D = Diagonal(d)
        mul!(Ms1, M, D)
        ldiv!(Ms2, D, Ms1)
        cond(Ms2)
    end
    n = size(M, 2)
    d0 = ones(real(eltype(M)), n)
    opts = Optim.Options(
        store_trace       = false,
        show_trace        = false,
        show_every        = 10,
        iterations        = 1000,
        allow_f_increases = false,
        time_limit        = 100,
        x_tol             = 0,
        f_abstol          = 0,
        g_tol             = tol,
        f_calls_limit     = 0,
        g_calls_limit     = 0,
    )

    res = Optim.optimize(
        scaleloss,
        d0,
        # i == 1 ? ParticleSwarm() : BFGS(alphaguess = LineSearches.InitialStatic(alpha=0.9), linesearch = LineSearches.HagerZhang()),
        n == 1 ? BFGS() : NelderMead(),
        opts,
    )
    d = res.minimizer
    abs.(d)
end

"""
    loop_scale(L::LTISystem, w = 0)

Find the optimal diagonal scaling matrix `D` such that `D\\L(iw)*D` has a minimized condition number at frequency `w`. Applicable to square `L` only.
Use [`loop_scaling`](@ref) to obtain `D`.
"""
function loop_scale(L::LTISystem, w = 0)
    M = freqresp(L, w)
    D = loop_scaling(M)
    ny = L.ny
    D = Diagonal(D)
    inv(D)*L*D
end

# function loop_scaling(M0::Matrix, tol=1e-4)
#     M = [
#         0I inv(M0)
#         M0     0I
#     ]
#     Ms1 = similar(M)
#     Ms2 = similar(M)
#     function scaleloss(d)
#         d02[1:n] .= d
#         d02[n+1:end] .= d
#         D = Diagonal(d02)
#         mul!(Ms1, M, D)
#         ldiv!(Ms2, D, Ms1)
#         cond(Ms2)
#     end
#     n = size(M0, 2)
#     d0 = ones(real(eltype(M)), n)
#     d02 = ones(real(eltype(M)), 2n)
#     opts = Optim.Options(
#         store_trace       = false,
#         show_trace        = false,
#         show_every        = 10,
#         iterations        = 1000,
#         allow_f_increases = false,
#         time_limit        = 100,
#         x_tol             = 0,
#         f_abstol          = 0,
#         g_tol             = tol,
#         f_calls_limit     = 0,
#         g_calls_limit     = 0,
#     )

#     res = Optim.optimize(
#         scaleloss,
#         d0,
#         # i == 1 ? ParticleSwarm() : BFGS(alphaguess = LineSearches.InitialStatic(alpha=0.9), linesearch = LineSearches.HagerZhang()),
#         n == 1 ? BFGS() : NelderMead(),
#         opts,
#     )
#     d = res.minimizer
#     abs.(d)
# end
