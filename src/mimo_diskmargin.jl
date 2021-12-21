using ControlSystems, IntervalArithmetic, MonteCarloMeasurements
using Optim

# TODO: introduce abstract type UncertaintyElement with subtypes interval uncertainty and ParticleUncertainty
# it's also beneficial to have real uncertainties (gain/parameter variation), phase uncertainties (abs = 1), and fully complex. Bounds are expected to get tighter if we use a more restrictive class of uncertainties.
# NOTE: for complex structured perturbations the typical upper bound algorithm appears to be quite tight. For real perturbations maybe the interval method has something to offer?
# NOTE: The naive algorithm using Optim appears to work quite well, but is a bit slow perhaps. The bounds computed using Interval methods do not seem to agree with the Optim method, in particular, they can be both above and below the Optim bound. Interestingly, the interval bound is never below the MCM bound. The MCM bound and the Optim bound are sometimes tight, but sometimes they differ a lot. So it could be that the interval bound is a true upper bound that happens to be better than the optim bound at times, but most often worse. Take the minimum of the two? The interval bound is a bit faster to compute, at least for small examples. Matlab's lower and upper bound on the test problem are tight, so the interval method must be wrong
# TODO: work in block structure in mussv
δ(a=1, N=0) = Complex(Interval(-a,a), Interval(-a,a))
function δp(a = 1, N=32)
    r = Particles(N, Uniform(0, a))
    θ = Particles(N, Uniform(0, 2pi))
    r*cis(θ)
end
Δ(n, a, δ = δ, N=32) = diagm([δ(a, N) for _ in 1:n])

function Base.:∈(a::Real, b::Complex{<:AbstractParticles})
    mi, ma = pextrema(b.re)
    mi <= a <= ma || return false
    mi, ma = pextrema(b.im)
    mi <= a <= ma
end


any0det(D::Matrix{<:Complex{<:Interval}}) = 0 ∈ det(D)

"""
    bisect_a(P, K, w; W = (:), Z = (:), au0 = 3.0, tol = 0.001, N = 32, upper = false)

For each frequency in `w`, find the largest a such that the loop with uncertainty elements of norm no greater than a, located at the inputs `W` and outputs `Z` of `P`, is stable.

By default, a conservative lower bound on the disk margin is returned. The conservatism comes from multiple factors
- Each uncertainty element is represented as an IntervalBox in the complex plane rather than a ball with radius `a`.
- IntervalArithmetic is inherently conservative due to the "dependency problem", i.e., `x-x ≠ 0`.

If `upper = true`, a Monte-Carlo approach is used to find an upper bound of the disk margin, this will be optimistic. If you calculate both and they are tight, you have a good indication of the true disk margin.

# Arguments:
- `P`: Plant
- `K`: Controller
- `w`: Frequency vector
- `W`: The inputs at which to insert a perturbation, defaults to all.
- `Z`: The outputs at which to insert a perturbation, defaults to all.
- `au0`: The largest a to try in the bisection. `a >= 2` yields infinite gain margin.
- `tol`: The tolerance in the bisection and the resolution of the resulting `a`.
- `N`: The number of samples for the upper bound estimation. 
- `upper`: Calculate upper or lower bound?
"""
function bisect_a(args...;  au0 = 3.0, tol=1e-3, kwargs...)
    M0, D = get_M(args...; kwargs...)
    iters = ceil(Int, log2(au0/tol))
    @views map(axes(M0, 3)) do i
        # @show i
        au = au0
        al = 0.0
        local a
        for j = 1:iters
            a = (au+al)/2
            if any0det(I - (a*M0[:,:,i])*D) # 0 ∈ det(I - (a*M0[:,:,i])*D)
                au = a
            else
                al = a
            end
        end
        al
    end
end

function get_M(P, K, w; W = (:), Z = (:), N = 32, upper=false)
    Z1 = W2 = Z == (:) ? (1:P.ny) : Z
    W1 = Z2 = W == (:) ? (1:P.nu) : W
    ny,nu = length(Z2), length(W2)
    if upper
        D = Δ(ny+nu, 1, δp, N)
    else
        D = Δ(ny+nu, 1, δ, N)
    end
    M = feedback(P, K; W2, Z2, Z1, W1)
    M0 = freqresp(M, w)
    M0 = permutedims(M0, (2,3,1))
    M0, D
end


# function muloss(M, d)
#     D = Diagonal(d)
#     opnorm(D\M*D)
# end


"""
    mussv(M; tol=1e-4)

Compute (an upper bound of) the structured singular value for diagonal Δ of complex perturbations.
`M` is assumed to be an (n × n × N_freq) array or a matrix.
"""
function mussv(M::AbstractArray{T}; tol=1e-4) where T
    Ms1 = similar(M[:,:,1])
    Ms2 = similar(M[:,:,1])
    function muloss(M, d)
        D = Diagonal(d)
        mul!(Ms1, M, D)
        ldiv!(Ms2, D, Ms1)
        opnorm(Ms2)
    end
    d0 = ones(real(T), size(M, 2))
    mu = map(axes(M, 3)) do i
        @views M0 = M[:,:,i]
        res = Optim.optimize(
            d->muloss(M0,d),
            d0,
            # i == 1 ? ParticleSwarm() : BFGS(alphaguess = LineSearches.InitialStatic(alpha=0.9), linesearch = LineSearches.HagerZhang()),
            i == 1 ? ParticleSwarm() : NelderMead(), # Initialize using Particle Swarm
            Optim.Options(
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
            ),
            # autodiff = :forward,
        )
        d0 = res.minimizer # update initial guess
        res.minimum
    end
    M isa AbstractMatrix ? mu[] : mu
end

"""
    sim_diskmargin(L, w::AbstractVector, σ::Real = 0)
    sim_diskmargin(L, σ::Real = 0)

Simultaneuous diskmargin at the outputs of `L`. 
"""
function sim_diskmargin(L,w::AbstractVector,σ::Real=0)
    # L = [ss(zeros(P.ny, P.ny)) P;-C ss(zeros(C.ny, C.ny))]
    # X = S+(σ-1)/2*I = lft([(1+σ)/2 -1;1 -1], L)
    n = L.ny
    X = ss(kron([(1+σ)/2 -1;1 -1], I(n)), L.timeevol)
    S̄ = feedback(X, L, U1 = n+1:2*n, Y1 = n+1:2*n, U2 = 1:n, Y2 = 1:n, Z1 = n+1:2*n, W1 = n+1:2*n, pos_feedback=true)
    M0 = permutedims(freqresp(S̄, w), (2,3,1))
    mu = mussv(M0)
    imu = inv.(mussv(M0))
    simultaneous = [Diskmargin(imu; ω0 = w, L) for (imu, w) in zip(imu,w)]
end

"""
    sim_diskmargin(P::LTISystem, C::LTISystem, σ::Real = 0)

Simultaneuous diskmargin at both outputs and inputs of `P`. 
"""
function sim_diskmargin(P::LTISystem, C::LTISystem, σ::Real=0)
    L = [ss(zeros(P.ny, P.ny)) P;-C ss(zeros(C.ny, C.ny))]
    sim_diskmargin(L,σ)
end

function sim_diskmargin(L,σ::Real=0)
    m = sim_diskmargin(L, LinRange(-3, 3, 500), σ)
    m = argmin(d->d.α, m)
end


"""
    ControlSystems.diskmargin(P::LTISystem, C::LTISystem, σ, w::AbstractVector, args...; kwargs...)

Simultaneuous diskmargin at outputs, inputs and input/output simultaneously of `P`. 
Returns a named tuple with the fields `input, output, simultaneous_input, simultaneous_output, simultaneous` where `input` and `output` represent loop-at-a-time margins, `simultaneous_input` is the margin for simultaneous perturbations on all inputs and `simultaneous` is the margin for perturbations on all inputs and outputs simultaneously.
"""
function ControlSystems.diskmargin(P::LTISystem, C::LTISystem, σ, w::AbstractVector, args...; kwargs...)
    L = C*P
    input = diskmargin(L, σ, w)
    simultaneous_input = sim_diskmargin(L,w,σ)
    L = P*C
    simultaneous_output = sim_diskmargin(L,w,σ)
    L = [ss(zeros(P.ny, P.ny)) P;-C ss(zeros(C.ny, C.ny))]
    simultaneous = sim_diskmargin(L,w,σ)
    (; input, output, simultaneous_input, simultaneous_output, simultaneous)
end

