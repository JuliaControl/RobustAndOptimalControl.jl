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


abstract type UncertainElement{C, F} end

struct δ{C, F} <: UncertainElement{C, F}
    val::C
    radius::F
    name::Symbol
    # function δr(val = 0, interval = Interval(val-1, val+1))
    #     lower ≤ val ≤ upper || throw(ArgumentError("val must be between lower and upper"))
    #     val, lower, upper = promote(val, lower, upper)
    #     new{typeof(val)}(val, lower, upper)
    # end
    function δ(val, radius, name = gensym())
        radius ≥ 0 || throw(ArgumentError("radius must be positive"))
        new{typeof(val), typeof(radius)}(val, radius, name)
    end
end

function δr(val::Real = 0, radius::Real = 1, args...)
    radius ≥ 0 || throw(ArgumentError("radius must be positive"))
    val, radius = promote(val, radius)
    δ(val, radius, args...)
end

function δc(val = 0, radius=1, args...)
    radius ≥ 0 || throw(ArgumentError("radius must be positive"))
    val = complex(val)
    δ(val, radius, args...)
end

Base.size(::δ) = (1,)
Base.size(::δ, n) = 1
Base.length(::δ) = 1
Base.eltype(::UncertainElement{C}) where C = C

function Base.promote_rule(::Type{N}, ::Type{δ{C, F}}) where {N <: Number, F, C}
    δ{promote_type(N, C), promote_type(N, F)}
end

Base.convert(::Type{δ{C, F}}, n::Number) where {C,F} = δ(n, 0)
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
    δ(*(n, d2.val), n*d2.radius, d2.name)
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
    # δr(d1.val+d2.val, d1.lower+d2.lower, d1.upper+d2.upper)
    uss(d1)*uss(d2)
end

function Base.:(/)(d1::δ, d2::δ)
    # δr(d1.val+d2.val, d1.lower+d2.lower, d1.upper+d2.upper)
    uss(d1)*inv(uss(d2))
end

function Base.:(/)(d1::Number, d2::δ)
    # δr(d1.val+d2.val, d1.lower+d2.lower, d1.upper+d2.upper)
    d1*inv(uss(d2))
end

function Base.inv(d2::δ)
    # δr(d1.val+d2.val, d1.lower+d2.lower, d1.upper+d2.upper)
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
# UncertainSS(sys::ExtendedStateSpace{TE}, Δ)

function Base.getproperty(sys::UncertainSS, s::Symbol)
    s ∈ fieldnames(typeof(sys)) && return getfield(sys, s)
    if s === :nz
        return foldl((l, d)->l + size(d, 1), sys.Δ, init=0)
    elseif s === :nw
        return foldl((l, d)->l + size(d, 2), sys.Δ, init=0)
    elseif s === :zinds
        return 1:sys.nz
    elseif s === :yinds
        return sys.nz .+ (1:size(sys.C2, 1))
    elseif s === :winds
        return 1:sys.nw
    elseif s === :uinds
        return sys.nw .+ (1:size(sys.B2, 2))
    end
    getproperty(getfield(sys, :sys), s)
end

function Base.promote_rule(::Type{U}, ::Type{δ{C, F}}) where {U <: UncertainSS, F, C}
    UncertainSS
end

function Base.promote_rule(::Type{U}, ::Type{L}) where {U <: UncertainSS, L <: LTISystem}
    UncertainSS
end

function Base.convert(::Type{U}, d::δ) where {U <: UncertainSS}
    uss(d)
end

function Base.convert(::Type{U}, s::LTISystem) where {U <: UncertainSS}
    s isa UncertainSS && return s
    sys = ss(s)
    # A,B,C,D = ssdata(sys)
    # esys = ss(A,[], B, [], C, [], [], [], D, sys.timeevol)
    esys = partition(sys, w=1:0, z = 1:0)
    UncertainSS(esys, [])
end

Base.:+(sys::ST, n::Union{δ, Number}) where ST <: UncertainSS = +(sys, uss(n))
Base.:+(n::Union{δ, Number}, sys::ST) where ST <: UncertainSS = +(uss(n), sys)

function uss(D11, D12, D21, D22, Δ, Ts=nothing)
    D11, D12, D21, D22 = vcat.((D11, D12, D21, D22))
    sys = ExtendedStateSpace(zeros(0,0), zeros(0,size(D21,2)), zeros(0,size(D12,2)), zeros(size(D12,1),0), zeros(size(D21,1),0), D11, D12, D21, D22, Ts)
    #                                       B1                   B2                    C1                        C2
    UncertainSS(sys, Δ)
end

function uss(d::δ{C,F}, Ts = nothing) where {C,F}
    # sys = partition(ss([d.val d.radius; 1 0]), u=2, y=2, w=1, z=1)
    uss(0, 1, d.radius, d.val, [d], Ts)
end

uss(n::Number, Ts=nothing) = uss(zeros(0,0), zeros(0,1), zeros(1,0), 1, [], Ts)

function uss(D::AbstractArray, Δ, Ts=nothing)
    length(Δ) == size(D,1) || throw(DimensionMismatch("length(Δ) != size(D,1)"))
    uss(D, zeros(size(D,1),0), zeros(0,size(D,2)), zeros(0,0), Δ, Ts)
end

uss(s::UncertainSS) = s

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
    if true
        sys = feedback(ss(sys1), ss(sys2),
            Z1=[s1.zinds; s1.yinds], Z2=s2.zinds,
            W2=[s2.winds; s2.uinds], W1=s1.winds,
            Y1=[], U2=[],
            U1=s1.uinds, Y2=s2.yinds,
            pos_feedback=true
        )


        sys = partition(sys, sys1.nz+sys2.nz, sys1.nw + sys2.nw)
        UncertainSS(sys, [s1.Δ; s2.Δ])
    else
        s1.nx == s2.nx == 0 || error("Not handled yet")
        M11, M12, M21, M22 = s1.D11, s1.D12, s1.D21, s1.D22
        Q11, Q12, Q21, Q22 = s2.D11, s2.D12, s2.D21, s2.D22
        A = [M11 M12*Q21; 0I Q11] # shall lower right be M11 or Q11? ss mult (Q11) differs from Zhou (M11)
        B = [M12*Q22; Q12]
        C = [M21 M22*Q21]
        D = M22*Q22
        uss(A,B,C,D, [s1.Δ; s2.Δ])
    end
end


function Base.:+(sys1::UncertainSS, sys2::UncertainSS) 
    s1 = sys1.sys
    s2 = sys2.sys
    s1.nx == s2.nx == 0 || error("Not handled yet")
    M11, M12, M21, M22 = s1.D11, s1.D12, s1.D21, s1.D22
    Q11, Q12, Q21, Q22 = s2.D11, s2.D12, s2.D21, s2.D22
    A = ControlSystems.blockdiag(M11, Q11)
    B = [M12; Q12]
    C = [M21 Q21]
    D = M22+Q22
    uss(A,B,C,D, [sys1.Δ; sys2.Δ])
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
        Δ = Diagonal(Δ)
    end

    if type === :l
        error("Invalid type of lft ($type), specify type=:u")
    else
        lft(ss(G.sys), ss(Δ), :u)
    end
end


## Uncertain dynamics

struct δDyn{F} <: UncertainElement{Complex{F}, F}
    ny::Int
    nu::Int
    bound::F
end

Base.size(d::δDyn, n=:) = (d.ny, d.nu)[n]
Base.length(d::δDyn) = d.ny + d.nu # NOTE: feels dangerous not to have prod(size)


function δss(ny, nu, Ts=nothing)
    # Δ = [δc(0, 1) for i = 1:(ny+nu)] # NOTE: perhaps not diagonal?
    Δ = δDyn(nu, ny, 1.0)
    D11 = Matrix([zeros(nu, ny) I(nu); I(ny) zeros(ny, nu)])
    # D22 = zeros(ny, nu)
    # D12 = zeros(ny+nu, nu)
    # D21 = zeros(ny, nu+ny)
    # uss(D11, D12, D21, D22, Δ, Ts)
    uss(D11, Δ, Ts)
end

## Sampling =======================================================
MonteCarloMeasurements.nominal(d::UncertainElement) = d.val

function Base.rand(d::δ{<:Real}, N::Int)
    Particles(N, Uniform(d.val - d.radius, d.val + d.radius))
end

function Base.rand(d::δ{<:Complex}, N::Int)
    r = Particles(N, Uniform(0, d.radius))
    θ = Particles(N, Uniform(0, 2pi))
    c = r*cis(θ)
    d.val == 0 ? c : c + d.val
end

normalize(d::δ) = δ(0*d.val, d.radius\d.radius, d.name)

function Base.rand(s::UncertainSS, N::Int)
    Δ = rand.(normalize.(s.Δ), N)
    lft(s, Δ, :u)
end

## MCM

function δr(N::Int)
    Particles(N, Uniform(-1, 1))
end

function δc(N::Int)
    r = Particles(N, Uniform(0, 1))
    θ = Particles(N, Uniform(0, 2pi))
    c = r*cis(θ)
end

Δ(n, δ) = diagm([δ() for _ in 1:n])
