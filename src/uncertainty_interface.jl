abstract type UncertainElement{C, F} end

struct δ{C, F<:Real} <: UncertainElement{C, F}
    val::C
    radius::F
    name::Symbol
    # function δr(val = 0, interval = Interval(val-1, val+1))
    #     lower ≤ val ≤ upper || throw(ArgumentError("val must be between lower and upper"))
    #     val, lower, upper = promote(val, lower, upper)
    #     new{typeof(val)}(val, lower, upper)
    # end
    function δ(val, radius, name = gensym("δ"))
        radius ≥ 0 || throw(ArgumentError("radius must be positive"))
        new{typeof(val), typeof(radius)}(val, radius, name)
    end
end

Base.zero(::Type{δ{C, F}}) where {C,F} = δ(zero(C), zero(F))

"""
    δr(val::Real = 0.0, radius::Real = 1.0, name)

Create a real, uncertain parameter. If no name is given, a boring name will be generated automatically.
"""
function δr(val::Real = 0.0, radius::Real = 1.0, args...)
    radius ≥ 0 || throw(ArgumentError("radius must be positive"))
    val, radius = promote(val, radius)
    δ(val, radius, args...)
end

"""
    δr(val::Real = 0.0, radius::Real = 1.0, name)

Create a complex, uncertain parameter. If no name is given, a boring name will be generated automatically.
"""
function δc(val = complex(0.0), radius=1.0, args...)
    radius ≥ 0 || throw(ArgumentError("radius must be positive"))
    val = complex(val)
    δ(val, radius, args...)
end

Base.size(::δ) = (1,)
Base.size(::δ, n) = 1
Base.length(::δ) = 1
Base.eltype(::UncertainElement{C}) where C = C
Base.eltype(::Type{<:UncertainElement{C}}) where C = C

function Base.promote_rule(::Type{N}, ::Type{δ{C, F}}) where {N <: Number, F, C}
    δ{promote_type(N, C), promote_type(N, F)}
end

Base.convert(::Type{δ{C, F}}, n::Number) where {C,F} = δ(n, 0.0)
Base.convert(::Type{δ{C, F}}, d::δ) where {C,F} = δ(C(d.val), F(d.radius), d.name)

#  + and - we define here, all other ops generate ExtendedStateSpace. THis is to avoid creating a complicated tagging mechanism
function Base.:(+)(d1::δ, d2::δ)
    if d1.name === d2.name
        δ(+(d1.val, d2.val), +(d1.radius, d2.radius), d1.name)
    else
        uss(d1) + uss(d2)
    end
end

function Base.:(+)(n::Number, d2::δ)
    δ(+(n, d2.val), d2.radius, d2.name)
end

Base.:(+)(d2::δ, n::Number) = +(n, d2)

function Base.:(*)(n::Number, d2::δ)
    δ(*(n, d2.val), abs(n)*d2.radius, d2.name)
end

function Base.:(-)(d2::δ)
    δ(-d2.val, d2.radius, d2.name)
end

function Base.:(-)(d1::δ, d2::δ)
    # δr(d1.val+d2.val, d1.lower+d2.lower, d1.upper+d2.upper)
    if d1.name === d2.name
        δ(-(d1.val, d2.val), -(d1.radius, d2.radius), d1.name)
    else
        # δ(-(d1.val, d2.val), +(d1.radius, d2.radius))
        uss(d1) + (-uss(d2))
    end
end

function Base.:(*)(d1::δ, d2::δ)
    uss(d1)*uss(d2)
end

function Base.:(/)(d1::δ, d2::δ)
    uss(d1)*inv(uss(d2))
end

function Base.:(/)(d1::Number, d2::δ)
    d1*inv(uss(d2))
end

function Base.inv(d2::δ)
    inv(uss(d2))
end

"""
    UncertainSS{TE} <: AbstractStateSpace{TE}

Represents LFT_u(M, Diagonal(Δ))
"""
struct UncertainSS{TE} <: AbstractStateSpace{TE}
    sys::ExtendedStateSpace{TE}
    Δ
end

function Base.getproperty(sys::UncertainSS, s::Symbol)
    s ∈ fieldnames(typeof(sys)) && return getfield(sys, s)
    if s === :nz
        # return foldl((l, d)->l + size(d, 1), sys.Δ, init=0)
        return size(sys.C1, 1)
    elseif s === :nw
        # return foldl((l, d)->l + size(d, 2), sys.Δ, init=0)
        return size(sys.B1, 2)
    elseif s === :zinds
        return 1:sys.nz
    elseif s === :yinds
        return sys.nz .+ (1:size(sys.C2, 1))
    elseif s === :winds
        return 1:sys.nw
    elseif s === :uinds
        return sys.nw .+ (1:size(sys.B2, 2))
    elseif s ===:M
        return sminreal(performance_mapping(sys.sys))
        # return sys.sys
    elseif s ===:delta
        return sys.Δ
    end
    getproperty(getfield(sys, :sys), s)
end

for f in [:system_mapping, :performance_mapping, :ss]
    @eval $f(s::UncertainSS, args...; kwargs...) = $f(s.sys, args...; kwargs...)
end

function Base.promote_rule(::Type{U}, ::Type{δ{C, F}}) where {U <: UncertainSS, F, C}
    UncertainSS
end

function Base.promote_rule(::Type{U}, ::Type{L}) where {U <: UncertainSS, L <: LTISystem}
    UncertainSS
end

function Base.promote_rule(::Type{U}, ::Type{L}) where {U <: UncertainSS, L <: AbstractArray}
    UncertainSS
end

function Base.convert(::Type{U}, d::δ) where {U <: UncertainSS}
    uss(d)
end

function Base.convert(::Type{U}, d::AbstractArray) where {U <: UncertainSS}
    uss(d)
end

function Base.convert(::Type{U}, s::LTISystem) where {U <: UncertainSS}
    s isa UncertainSS && return s
    sys = ss(s)
    esys = partition(sys, 0, 0)
    UncertainSS(esys, [])
end

Base.:+(sys::UncertainSS{TE}, n::δ) where TE <: ControlSystems.TimeEvolution = +(sys, uss(n))
Base.:+(n::δ, sys::UncertainSS{TE}) where TE <: ControlSystems.TimeEvolution = +(uss(n), sys)

Base.:+(sys::UncertainSS{TE}, n::Number) where TE <: ControlSystems.TimeEvolution = +(sys, uss(n))
Base.:+(n::Number, sys::UncertainSS{TE}) where TE <: ControlSystems.TimeEvolution = +(uss(n), sys)

Base.:*(G::LTISystem, d::UncertainElement) = uss(ss(G)) * uss(d)
Base.:*(d::UncertainElement, G::LTISystem) = uss(d) * uss(ss(G))

"""
    uss(D11, D12, D21, D22, Δ, Ts = nothing)

Create an uncertain statespace object with only gin matrices.
"""
function uss(D11, D12, D21, D22, Δ, Ts=nothing)
    D11, D12, D21, D22 = vcat.((D11, D12, D21, D22))
    sys = ExtendedStateSpace(zeros(0,0), zeros(0,size(D21,2)), zeros(0,size(D12,2)), zeros(size(D12,1),0), zeros(size(D21,1),0), D11, D12, D21, D22, Ts)
    #                                       B1                   B2                    C1                        C2
    UncertainSS(sys, Δ)
end

"""
    uss(d::δ{C, F}, Ts = nothing)

Convert a δ object to an UncertainSS
"""
function uss(d::δ{C,F}, Ts = nothing) where {C,F}
    # sys = partition(ss([d.val d.radius; 1 0]), u=2, y=2, w=1, z=1)
    uss(0, 1, d.radius, d.val, [normalize(d)], Ts)
end

"""
    uss(d::AbstractVector{<:δ}, Ts = nothing)

Create a diagonal uncertain statespace object with the uncertain elements `d` on the diagonal.
"""
function uss(d::AbstractVector{<:δ}, Ts = nothing)
    vals = getfield.(d, :val)
    rads = getfield.(d, :radius)
    n = length(d)
    uss(0I(n), I(n), diagm(rads), diagm(vals), normalize.(d), Ts)
end

uss(n::Number, Ts=nothing) = uss(zeros(0,0), zeros(0,1), zeros(1,0), 1, [], Ts)

"""
    uss(D::AbstractArray, Δ, Ts = nothing)

If only a single `D` matrix is provided, it's treated as `D11` if Δ is given, and as `D22` if no Δ is provided.
"""
function uss(D::AbstractArray, Δ=[], Ts=nothing)
    if length(Δ) == size(D,1)
        uss(D, zeros(size(D,1),0), zeros(0,size(D,2)), zeros(0,0), Δ, Ts)
    elseif isempty(Δ)
        uss(zeros(0,0), zeros(0,size(D,2)), zeros(size(D,1),0), D, Δ, Ts)
    else
        throw(DimensionMismatch("length(Δ) != size(D,1)"))
    end
end

uss(s::UncertainSS) = s
uss(s::AbstractStateSpace) = convert(UncertainSS, s)

# Lower LFT
# function Base.inv(s::UncertainSS)
#     s.nx == 0 || error("Haven't considered this case yet")
#     @unpack D11, D12, D21, D22 = s
#     iD11 = inv(D11)
#     sys = uss(iD11, -iD11*D12, D21*iD11, D22-D21*iD11*D12, s.Δ, s.timeevol)
# end

function Base.inv(s::UncertainSS) # Upper LFT
    s.nx == 0 || error("Haven't considered this case yet")
    @unpack D11, D12, D21, D22 = s
    iD22 = inv(D22)
    sys = uss(D11-D12*iD22*D21, D12*iD22, -iD22*D21, iD22, s.Δ, s.timeevol)
end



function Base.:*(s1::UncertainSS, s2::UncertainSS) 
    sys1 = s1.sys
    sys2 = s2.sys
    sys = invert_mappings(invert_mappings(sys1)*invert_mappings(sys2))
    UncertainSS(sys, [s1.Δ; s2.Δ])
end

function Base.:*(s1::UncertainSS, n::Number)  # reverse case is handled by switching
    sys = invert_mappings(n*invert_mappings(s1.sys))
    UncertainSS(sys, s1.Δ)
end


function Base.:+(s1::UncertainSS, s2::UncertainSS) 
    sys1 = s1.sys
    sys2 = s2.sys
    sys = invert_mappings(invert_mappings(sys1)+invert_mappings(sys2))
    UncertainSS(sys, [s1.Δ; s2.Δ])
end

function Base.:-(s::UncertainSS) 
    UncertainSS(-s.sys, s.Δ)
end

function Base.hcat(systems::UncertainSS...)
    # Perform checks
    timeevol = common_timeevol(systems...)

    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same ouput dimension")
    end

    A = blockdiag([s.A for s in systems]...)

    B1 = blockdiag([s.B1 for s in systems]...)
    B2 = blockdiag([s.B2 for s in systems]...)

    C2 = reduce(hcat, [s.C2 for s in systems])
    C1 = blockdiag([s.C1 for s in systems]...)

    D21 = reduce(hcat, [s.D21 for s in systems])
    D22 = reduce(hcat, [s.D22 for s in systems])
    D11 = blockdiag([s.D11 for s in systems]...)
    D12 = blockdiag([s.D12 for s in systems]...)

    sysnew = ss(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
    return UncertainSS(sysnew, reduce(vcat, [s.Δ for s in systems]))
end

function Base.vcat(systems::UncertainSS...)
    # Perform checks
    timeevol = common_timeevol(systems...)

    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end

    A = blockdiag([s.A for s in systems]...)

    B2 = reduce(vcat, [s.B2 for s in systems])
    B1 = blockdiag([s.B1 for s in systems]...)

    C1 = blockdiag([s.C1 for s in systems]...)
    C2 = blockdiag([s.C2 for s in systems]...)

    D12 = reduce(vcat, [s.D12 for s in systems])
    D22 = reduce(vcat, [s.D22 for s in systems])
    D11 = blockdiag([s.D11 for s in systems]...)
    D21 = blockdiag([s.D21 for s in systems]...)

    sysnew = ss(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
    return UncertainSS(sysnew, reduce(vcat, [s.Δ for s in systems]))
end

function ControlSystems.lft(G::UncertainSS, Δ::AbstractArray=G.Δ, type=:u)
    if ndims(Δ) == 1 
        Δ = append(ss.(Δ)...)
    end

    if type !== :u
        error("Invalid type of lft ($type), specify type=:u")
    else
        lft(ss(G.sys), ss(Δ), :u)
    end
end

function ControlSystems.lft(G::UncertainSS, K::LTISystem, type=:l)
    if type !== :l
        error("Invalid type of lft ($type), specify type=:l")
    else
        sys = lft(ss(G.sys), ss(K), :l)
        sys = partition(sys, sys.nu, sys.nu)
        UncertainSS(sys, G.Δ)
    end
end

## Uncertain dynamics

struct δDyn{F, TE} <: UncertainElement{Complex{F}, F}
    ny::Int
    nu::Int
    nx::Int
    bound::F
    timeevol::TE
    name::Symbol
end

function Base.getproperty(d::δDyn, s::Symbol)
    s ∈ fieldnames(typeof(d)) && return getfield(d, s)
    if s === :Ts
        return d.timeevol.Ts
    else
        throw(ArgumentError("$(typeof(d)) has no property named $s"))
    end
end


function δDyn(ny::Int, nu::Int; bound=1, nx::Int = 2, Ts=0, name=gensym("δDyn"))
    timeevol = Ts > 0 ? Discrete(Ts) : Continuous()
    δDyn(ny,nu,nx,bound,timeevol,name)
end

Base.size(d::δDyn, n=:) = (d.ny, d.nu)[n]
# Base.length(d::δDyn) = d.ny + d.nu # NOTE: feels dangerous not to have prod(size)

function δss(ny, nu; bound=1, kwargs...)
    Δ = δDyn(ny, nu; bound, kwargs...)
    b = sqrt(bound)
    # D11 = Matrix([zeros(nu, ny) b*I(nu); b*I(ny) zeros(ny, nu)])
    # D22 = zeros(ny,nu)
    # D12 = zeros(ny+nu, nu)
    # D21 = zeros(ny, nu+ny)
    D11 = zeros(nu, ny)
    D22 = zeros(ny, nu)
    D12 = b*I(nu)
    D21 = b*I(ny)
    uss(D11, D12, D21, D22, [normalize(Δ)], Δ.timeevol)
end

## Sampling =======================================================
MonteCarloMeasurements.nominal(d::UncertainElement) = d.val

function Base.rand(d::δ{V}, N::Int) where V <: Real
    T = float(V)
    d.radius == 0 && return Particles{T, N}(fill(T(d.val), N))
    Particles(N, Uniform(d.val - d.radius, d.val + d.radius))
end

function Base.rand(D::AbstractArray{<:δ}, N::Int)
    M = rand.(D, N)
    T = eltype(M[1])
    T.(M)
end

function Base.rand(d::δ{Complex{F}}, N::Int) where F
    T = float(F)
    d.radius == 0 && return Particles(fill(T(d.val), N))
    r = Particles(N, Uniform(0, d.radius))
    θ = Particles(N, Uniform(0, 2pi))
    c = r*cis(θ)
    d.val == 0 ? c : c + d.val
end

function Base.rand(d::δDyn, N::Int)
    G = map(1:N) do i
        while true
            G = ssrand(d.ny, d.nu, d.nx; Ts = d.timeevol isa Continuous ? nothing : d.Ts)
            try
                n = norm(G, Inf)
                return (1/n) * G
            catch
            end
        end
    end
    ss2particles(G)
end

normalize(d::δ) = δ(0*d.val, one(d.radius), d.name)
normalize(d::δDyn) = δDyn(d.ny, d.nu, d.nx, 1, d.timeevol, d.name)


function Base.rand(s::UncertainSS, N::Int)
    Δ = rand.(normalize.(s.Δ), N)
    lft(s, Δ, :u)
end

## MCM

# function LinearAlgebra.det(D::AbstractMatrix{<:Complex{<:AbstractParticles}})
#     D0 = similar(D, ComplexF64)
#     parts = map(1:nparticles(D[1].re)) do i
#         for j in eachindex(D0)
#             D0[j] = Complex(D[j].re.particles[i], D[j].im.particles[i])
#         end
#         det(D0)
#     end
#     Complex(StaticParticles(getfield.(parts, :re)), StaticParticles(getfield.(parts, :im)))
# end

function ControlSystems.tzeros(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractMatrix{T}, D::AbstractMatrix{T}) where T <: AbstractParticles
    bymap(tzeros, A, B, C, D)
end


using MonteCarloMeasurements: vecindex
function sys_from_particles(P, i)
    A,B,C,D = ssdata(P)
    ss(vecindex(A, i), vecindex(B, i), vecindex(C, i), vecindex(D, i))
end


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


function ss2particles(G::Vector{<:AbstractStateSpace})
    pdp(x) = Particles(permutedims(x))
    A = reduce(hcat, vec.(getproperty.(G, :A))) |> pdp
    B = reduce(hcat, vec.(getproperty.(G, :B))) |> pdp
    C = reduce(hcat, vec.(getproperty.(G, :C))) |> pdp
    D = reduce(hcat, vec.(getproperty.(G, :D))) |> pdp
    (; nx,ny,nu) = G[1]
    A = reshape(A, nx, nx)
    B = reshape(B, nx, nu)
    C = reshape(C, ny, nx)
    D = reshape(D, ny, nu)
    ss(A,B,C,D, G[1].timeevol)
end

function δr(N::Int)
    Particles(N, Uniform(-1, 1))
end

function δc(N::Int)
    r = Particles(N, Uniform(0, 1))
    θ = Particles(N, Uniform(0, 2pi))
    c = r*cis(θ)
end

Δ(n, δ) = Diagonal([δ() for _ in 1:n])

## Intervals
IntervalArithmetic.Interval(d::δ{R}) where R <: Real = 
    Interval(d.val-d.radius, d.val+d.radius)

function IntervalArithmetic.Interval(d::δ{C}) where C <: Complex 
    re = Interval(d.val.re-d.radius, d.val.re+d.radius)
    im = Interval(d.val.im-d.radius, d.val.im+d.radius)
    Complex(re, im)
end



"""
    block_structure(Δ)

Take a vector of uncertain elements and return a vector of vectors with the block
structure of perturbation blocks as described expected my μ-tools, i.e.
- `[-N, 0]` denots a repeated real block of size `N`
- `[N, 0]` denots a repeated complex block of size `N`
- `[ny, nu]` denots a full complex block of size `ny × nu`
"""
function block_structure(D)
    perm = sortperm(D, by=d->d.name)
    D = D[perm]
    blocks = Vector{Int}[]
    blockstart = 1
    for i in 2:length(D)
        if D[i].name == D[i-1].name
            if blockstart == 0
                blockstart = i-1
            end
        else
            push!(blocks, makeblock(D[blockstart:i-1]))
            blockstart = i
        end
    end
    push!(blocks, makeblock(D[blockstart:length(D)]))
    blocks, perm
end

function makeblock(d) # TODO: full uncertain block that is not a dynamical system is not supported
    N = length(d)
    if N > 1 || (N==1 && d[] isa δ) # Δ = δI -> μ(M) = ρ(M) 8.80
        if eltype(d[1]) <: Complex
            b = [N, 0] # complex*I
        else
            b = [-N, 0] # real*I
        end
    else
        b = [size(d[], 1), size(d[], 2)]
    end
end


function blocksort(P::UncertainSS)
    blocks, perm = block_structure(P.Δ)
    Poutinds = []
    Pininds = []
    # Build indices into M that corresponds to the Δ blocks and permute those indices with the same permutation as was applied to Δ
    for (i, d) in enumerate(P.Δ)
        if i == 1
            push!(Poutinds, 1:size(d, 2))
            push!(Pininds, 1:size(d, 1))
        else
            push!(Poutinds, (1:size(d, 2)) .+ Poutinds[i-1][end])
            push!(Pininds, (1:size(d, 1)) .+ Pininds[i-1][end])
        end
    end
    Poutinds = reduce(vcat, Poutinds[perm])
    Pininds = reduce(vcat, Pininds[perm])
    M = P.M[Poutinds, Pininds]
    blocks, M
end

function getinds(blocks::Vector{Vector{Int}}, d::Int)
    lengths = [b[2] <= 1 ? 1 : b[d] for b in blocks]
    endinds = cumsum(lengths)
    startinds = [1; endinds[1:end-1] .+ 1]
    inds = UnitRange.(0 .+ startinds, lengths .+ startinds .- 1)
end

# function getinds(D::AbstractVector, d::Int)
#     lengths = size.(D, d)
#     endinds = cumsum(lengths)
#     startinds = [1; endinds[1:end-1] .+ 1]
#     inds = UnitRange.(0 .+ startinds, lengths .+ startinds .- 1)
# end

# struct BlockDiag{T}
#     D::Vector{T}
#     rinds::Vector{UnitRange{Int}}
#     cinds::Vector{UnitRange{Int}}
# end
# BlockDiag(D::AbstractVector) = BlockDiag(D, getinds(D, 1), getinds(D, 2))

# Base.size(B::BlockDiag) = (B.rinds[end][end], B.cinds[end][end])

# function Base.:*(A, B::BlockDiag)
#     rinds = B.rinds
#     submats = map(enumerate(rinds)) do (i,inds)
#         A[:, inds] * B.D[i]
#     end
#     reduce(hcat, submats)
# end

# function quadform_minus(Q::BlockDiag, M, b, op)
#     # con = M'Q*M - b^2*Q ⪯ I(n)
#     rinds = Q.rinds
#     cinds = Q.cinds
#     map(eachindex(rinds)) do i
#         r = rinds[i]
#         c = cinds[i]
#         q = Q.D[i]
#         op(M[r,:]'q*M[c,:] .- b^2*q, I(size(M,2)))
#     end
# end



#=
D_from_Delta: I-full, full-I, diag-diag
Version that constructs from initial guess
=#


# See Skogestad 8.8.3
# for complex δ:
# Δ = δI -> μ(M) = ρ(M) 8.80
# Full block complex -> μ(M) = σ̄(M) = opnorm(M)
# In general ρ(M) <= μ(M) <= σ̄(M)

# function make_D(blocks)
#     D = []
#     for b in blocks
#         if b[2] == 0
#             n = abs(b[1])
#             T = b[1] > 0 ? ComplexF64 : Float64
#             push!(D, zeros(T, n, n))
#         else
#             push!(D, 0*I)
#         end
#     end
#     BlockDiagonal(D)
# end