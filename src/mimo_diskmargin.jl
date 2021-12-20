using ControlSystems, IntervalArithmetic, MonteCarloMeasurements

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

function LinearAlgebra.det(D::AbstractMatrix{<:Complex{<:AbstractParticles}})
    D0 = similar(D, ComplexF64)
    parts = map(1:nparticles(D[1].re)) do i
        for j in eachindex(D0)
            D0[j] = Complex(D[j].re.particles[i], D[j].im.particles[i])
        end
        det(D0)
    end
    Complex(StaticParticles(getfield.(parts, :re)), StaticParticles(getfield.(parts, :im)))
end
# M = [0 1; -0.1 -0.1]
# D = Δ(2, 1)
# 0 ∈ det(I-M*D)
# This computation must be done for all frequncies, since M = M(ω), so, for each frequency, bisect over a to find μ(ω) this sounds very expensive?

# 1. Construct MΔ system
# sys1 is both P and K, while sys2 is the deltas, ny+nu
# z1 are the outputs of P and w1 are the inputs of K
# if we set delta t0 0, we should get back the nominal closed-loop system
# P = ssrand(2,3,2, proper=true)
# K = ssrand(3,2,2, proper=true)


any0det(D::Matrix{<:Complex{<:Interval}}) = 0 ∈ det(D)
function any0det(D::Matrix{<:Complex{<:AbstractParticles}})
    D0 = similar(D, ComplexF64)
    maxre = maxim = -1
    minre = minim = 1
    for i = 1:nparticles(D[1].re)
        for j in eachindex(D0)
            D0[j] = Complex(D[j].re.particles[i], D[j].im.particles[i])
        end
        d = det(D0)
        maxre = max(maxre, d.re)
        minre = min(minre, d.re)
        maxim = max(maxim, d.im)
        minim = min(minim, d.im)
        if maxre > 0 && minre < 0 && maxim > 0 && minim < 0
            return true
        end
    end
    false
end

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


##

a = 10
P = [
        tf([1,-a^2], [1, 0, a^2]) tf([a, a], [1, 0, a^2])
        -tf([a, a], [1, 0, a^2]) tf([1,-a^2], [1, 0, a^2])
    ]
P = minreal(ss(P))
K = ss(1.0I(2))


ny,nu = size(P)
sys2 = ss(I(ny+nu)) # this formulation makes sense if sys2 is I + a*δ 

# sys1 = ControlSystems.append(P,K)
# z1 = 1:ny
# w1 = (1:ny) .+ ny
# M = feedback(sys1, sys2; Z1 = z1, W1 = w1)

# M = feedback(P, K, W2 = :, Z2 = :) # output everything
# feedback(M, ss(0*I(ny+nu)))

# 1. form system that exposes all inputs and outputs but also has feedback
M = feedback(P, K, W2 = :, Z2 = :)
@test poles(M) ≈ poles(feedback(P*K))
@test size(M) == (ny+nu, ny+nu)
# @test minreal(feedback(M, 0*sys2)) ≈ M


w = 2π .* exp10.(LinRange(-2, 2, 300))
# @time bisect_a(P, K, w)
##

# break at input (pass outputs through)
a = bisect_a(P, K, w; Z = [], tol=1e-4)
@test minimum(a) < 0.0998

# break at output (pass inputs through)
a = bisect_a(P, K, w; W = [], tol=1e-4)
@test minimum(a) < 0.0998

# break at both input and output
a = bisect_a(P, K, w; tol=1e-4)
au = bisect_a(P, K, w; tol=1e-4, N=640, upper=true)
@test minimum(a) < 0.0499


plot(w, a, xscale=:log10, xlabel="Frequency", ylims=(0,3))
plot!(w, au, xscale=:log10, xlabel="Frequency", ylims=(0,3))


##
w = 2π .* exp10.(LinRange(-2, 2, 300))


function vec2sys(v, ts = nothing)
    nx = round(Int, v[1])
    nu = round(Int, v[2])
    ny = round(Int, v[3])
    ai = (1:nx^2) .+ 3
    bi = (1:nx*nu) .+ ai[end]
    ci = (1:nx*ny) .+ bi[end]
    di = (1:nu*ny) .+ ci[end]
    A = reshape(v[ai], nx, nx)
    B = reshape(v[bi], nx, nu)
    C = reshape(v[ci], ny, nx)
    D = reshape(v[di], ny, nu)
    ts === nothing ? ss(A, B, C, D) : ss(A, B, C, D, ts)
end

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end

a = bisect_a(L3, ss(I(3), L3.Ts), w; tol=2e-3)
au = bisect_a(L3, ss(I(3), L3.Ts), w; tol=2e-3, upper=true, N=256)
plot(w, a, xscale=:log10, xlabel="Frequency", ylims=(0,Inf))
plot!(w, au, xscale=:log10, xlabel="Frequency", ylims=(0,Inf))


dm = diskmargin(L3, 0, w)
plot(w, dm[1,:])


dm = diskmargin(L3, 0, 4.05)
@test dm[1].α ≈ 0.794418036911981 rtol=1e-3


dm = diskmargin(L3, ss(I(3), L3.Ts), 0, w)


## Test diskmargin with particles in the system

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]*(1 + 0.1*Particles(32))
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


unsafe_comparisons(true)
dm = diskmargin(L3, 0, 4.05)



##

using Optim, Optim.LineSearches

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


P = L3
K = ss(I(3), L3.timeevol)
w = 2π .* exp10.(LinRange(-2, 2, 300))
# break at input (pass outputs through)
M,_ = get_M(P, K, w; Z = :)
M0 = M[:,:,100]

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

function get_worst_perturbation(M::AbstractMatrix{T}; tol=1e-4) where T
    Ms1 = similar(M)
    Ms2 = similar(M)
    function muloss(d)
        D = Diagonal(d)
        mul!(Ms1, M, D)
        ldiv!(Ms2, D, Ms1)
        opnorm(Ms2)
    end
    d0 = ones(real(T), size(M, 2))
    res = Optim.optimize(
        d->muloss(d),
        d0,
        ParticleSwarm(),
        Optim.Options(
            store_trace       = false,
            show_trace        = false,
            show_every        = 10,
            iterations        = 1500,
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
    res.minimizer
    # This needs to be transformed back if L is simultaneous on both inputs and outputs
end

@time mu = [mussv(M) for i = 1:10]
# mum = minimum(mu)
plot(w, mu, xscale=:log10)

a = bisect_a(P, K, w; Z = [], tol=1e-4)
au = bisect_a(P, K, w; Z = [], tol=1e-4, upper=true, N=2560)
plot!(w, inv.(a), xscale=:log10)
plot!(w, inv.(au), xscale=:log10)


##
function sim_diskmargin(L,w::AbstractVector,σ::Real=0)
    # L = [ss(zeros(P.ny, P.ny)) P;-C ss(zeros(C.ny, C.ny))]
    # X = S+(σ-1)/2*I = lft([(1+σ)/2 -1;1 -1], L)
    n = L.ny
    X = ss(kron([(1+σ)/2 -1;1 -1], I(n)), L.timeevol)
    S̄ = feedback(X, L, U1 = n+1:2*n, Y1 = n+1:2*n, U2 = 1:n, Y2 = 1:n, Z1 = n+1:2*n, W1 = n+1:2*n, pos_feedback=true)
    M0 = permutedims(freqresp(S̄, w), (2,3,1))
    mu = mussv(M0)
    imu = inv.(mussv(M0))
    @show get_worst_perturbation(M0[:,:,argmax(mu)])
    simultaneous = [Diskmargin(imu; ω0 = w, L) for (imu, w) in zip(imu,w)]
end

function sim_diskmargin(P::LTISystem, C::LTISystem, σ::Real=0)
    L = [ss(zeros(P.ny, P.ny)) P;-C ss(zeros(C.ny, C.ny))]
    sim_diskmargin(L,σ)
end

function sim_diskmargin(L,σ::Real=0)
    m = sim_diskmargin(L, LinRange(-3, 3, 500), σ)
    m = argmin(d->d.α, m)

end


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
C = K
dm = diskmargin(L3, 1.0*K, 0, w;)
plot(dm.simultaneous) # TODO: verkar ge fel svar, ungefär samma form som matlab men inte rätt värden

plot(dm.input)



## 
# NOTE: SISO och loop at a time blir rätt, men inte simultaneous. mussv verkar rätt, så kan vara fel på get_M
a = [-0.2 10;-10 -0.2]; b = I(2); c = [1 8;-10 1];
P = ss(a,b,c,0);
K = ss([1 -2;0 1]);
dm = diskmargin(K*P) # disk margins at plant inputs
dm = diskmargin(P*K); # disk margins at plant outputs
MMIO = diskmargin(P,K,w)

plot(MMIO.simultaneous)

w = 2π .* exp10.(LinRange(-3, 3, 300))

##
MMO = diskmargin(P,K,0,w)

plot(MMO.simultaneous, lab="simultaneous")
plot!(MMO.simultaneous_output, lab="output") # denna är rätt för små frekvenser men 2x fel för höga
plot!(MMO.simultaneous_input, lab="input") # samma för denna




L = K*P
M0 = permutedims(freqresp(feedback(L), w), (2,3,1))
mu = mussv(M0)
imu = inv.(mussv(M0))
simultaneous = [Diskmargin(imu; ω0 = w, L) for (imu, w) in zip(imu,w)]

plot(w, mu)

using MATLAB
mum = mat"mussv($M0, ones(2,2))"
plot(mu)
plot!(mum[1,:,:]')


