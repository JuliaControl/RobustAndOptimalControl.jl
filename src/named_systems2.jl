import ControlSystems as CS
import ControlSystems: nstates, blockdiag

struct NamedStateSpace{T,S} <: AbstractStateSpace{T} where S <: AbstractStateSpace{T}
    sys::S
    x
    u
    y
end

function Base.promote_rule(::Type{U}, ::Type{NamedStateSpace{T, S}}) where
    {T, U<:AbstractStateSpace{T} , S<:AbstractStateSpace{T}} 
    @show U,S
    @show inner = promote_type(U,S)
    NamedStateSpace{T, inner}
end

function Base.convert(::Type{NamedStateSpace{T, S}}, s::S) where {T, S <: AbstractStateSpace}
    named_ss(s, x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

function Base.convert(::Type{NamedStateSpace{T, S}}, s::U) where {T, S <: AbstractStateSpace, U <: AbstractStateSpace}
    s2 = Base.convert(S, s)
    named_ss(s2, x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

function Base.convert(::Type{NamedStateSpace{T, S}}, s::NamedStateSpace{T, U}) where {T, S <: AbstractStateSpace, U <: AbstractStateSpace}
    sys = Base.convert(S, s.sys)
    NamedStateSpace{T,typeof(sys)}(sys, s.x, s.u, s.y)
end

# Base.convert(::Type{RobustAndOptimalControl.NamedStateSpace{T, S}}, s::S) where {T, S<:RobustAndOptimalControl.NamedStateSpace} = s



const NamedIndex = Union{Symbol, Vector{Symbol}}

function Base.getproperty(G::NamedStateSpace, s::Symbol)
    s ∈ fieldnames(NamedStateSpace) && (return getfield(G,s))
    return getproperty(G.sys, s)
end

ControlSystems.numeric_type(G::NamedStateSpace) = ControlSystems.numeric_type(G.sys)

maybe_expand(s::Symbol, n::Int) = n == 1 ? [s] : [Symbol(string(s)*string(i)) for i in 1:n]
maybe_expand(v, n::Int) = v

function named_ss(sys::AbstractStateSpace{T};
    x = [Symbol("x$i") for i in 1:sys.nx],
    u = [Symbol("u$i") for i in 1:sys.nu],
    y = [Symbol("y$i") for i in 1:sys.ny],
    ) where T
    x = maybe_expand(x, sys.nx)
    u = maybe_expand(u, sys.nu)
    y = maybe_expand(y, sys.ny)
    length(x) == sys.nx ||
        throw(ArgumentError("Length of state names must match sys.nx ($(sys.nx))"))
    length(u) == sys.nu ||
        throw(ArgumentError("Length of input names must match sys.nu ($(sys.nu))"))
    length(y) == sys.ny ||
        throw(ArgumentError("Length of output names must match sys.ny ($(sys.ny))"))

    NamedStateSpace{T, typeof(sys)}(sys, x, u, y)
end

macro check_unique(ex, s = string(ex))
    quote
        $(esc(__source__))
        u = unique($(esc(ex)))
        if length(u) != length($(esc(ex)))
            rep = setdiff($(esc(ex)), u)
            throw(ArgumentError($(s)*" not unique. Repeated names: "*string(u)))
        end
    end
end

macro check_all_unique(s1, s2)
    quote
        vals = [getproperty($(esc(s1)), :x); getproperty($(esc(s2)), :x)]
        @check_unique vals "x"
        vals = [getproperty($(esc(s1)), :u); getproperty($(esc(s2)), :u)]
        @check_unique vals "u"
        vals = [getproperty($(esc(s1)), :y); getproperty($(esc(s2)), :y)]
        @check_unique vals "y"
    end
end

iterable(s::Symbol) = [s]
iterable(v) = v

function Base.getindex(sys::NamedStateSpace{T,S}, i::NamedIndex, j::NamedIndex) where {T,S}
    i,j = iterable.((i, j))
    ii = findall(sys.y .∈ (i, )) # findall(i .∈ (sys.y, ))
    jj = findall(sys.u .∈ (j, )) # findall(j .∈ (sys.u, ))

    return NamedStateSpace{T,S}(
        sys.sys[ii, jj],
        sys.x,
        sys.u[jj],
        sys.y[ii],
    )
end

function Base.getindex(sys::NamedStateSpace{T,S}, inds...) where {T,S}
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = CS.index2range(inds...)
    return NamedStateSpace{T,S}(
        sys.sys[rows, cols],
        sys.x,
        sys.u[cols],
        sys.y[rows],
    )
end

function Base.show(io::IO, G::NamedStateSpace)
    print(io, "Named")
    show(io, G.sys)
    print(io, "\nWith state  names: "); println(io, join(G.x, ' '))
    print(io, "     input  names: "); println(io, join(G.u, ' '))
    print(io, "     output names: "); println(io, join(G.y, ' '))
end

function Base.:-(s1::NamedStateSpace{T,S}) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        -s1.sys,
        s1.x,
        s1.u,
        s1.y,
    )
end

function Base.:+(s1::NamedStateSpace{T,S}, s2::NamedStateSpace{T,S}) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        s1.sys+s2.sys,
        [s1.x; s2.x],
        s1.u,
        s1.y,
    )
end

function Base.:*(s1::NamedStateSpace{T}, s2::NamedStateSpace{T}) where {T <: CS.TimeEvolution}
    @check_all_unique s1 s2
    if s1.u != s2.y
        connection_map = join(["$y -> $u" for (u,y) in zip(s1.u, s2.y) if u != y], '\n')
        @warn "Connected signals have different names\n $connection_map"
    end
    sys = s1.sys*s2.sys
    S = typeof(sys)
    return NamedStateSpace{T,S}(
        sys,
        [s1.x; s2.x],
        s2.u,
        s1.y,
    )
end
##

function Base.hcat(systems::NamedStateSpace{T,S}...) where {T,S}
    x = reduce(vcat, getproperty.(systems, :x))
    u = reduce(vcat, getproperty.(systems, :u))
    @check_unique x
    @check_unique u
    return NamedStateSpace{T,S}(
        hcat(getproperty.(systems, :sys)...),
        x,
        u,
        systems[1].y,
    )
end

function Base.vcat(systems::NamedStateSpace{T,S}...) where {T,S}
    x = reduce(vcat, getproperty.(systems, :x))
    y = reduce(vcat, getproperty.(systems, :y))
    @check_unique x
    @check_unique y
    return NamedStateSpace{T,S}(
        vcat(getproperty.(systems, :sys)...),
        x,
        systems[1].u,
        y,
    )
end


function ControlSystems.feedback(s1::NamedStateSpace{T,S}, s2::NamedStateSpace{T,S}; r = s1.u) where {T <: CS.TimeEvolution, S}
    sys = feedback(s1.sys, s2.sys)
    r = iterable(r)
    x1  = [s1.x; s2.x]
    @check_unique x1
    return NamedStateSpace{T,S}(sys, x1, r, s1.y)
end

function ExtendedStateSpace(P::NamedStateSpace; z=[], y=[], w=[], u=[])
    zi = [findfirst(==(zi), P.y) for zi in z]
    yi = [findfirst(==(yi), P.y) for yi in y]
    wi = [findfirst(==(wi), P.u) for wi in w]
    ui = [findfirst(==(ui), P.u) for ui in u]
    ss(P.A, P.B[:, wi], P.B[:, ui], P.C[zi, :], P.C[yi, :], 
        P.D[zi, wi], P.D[zi, ui], P.D[yi, wi], P.D[yi, ui], P.timeevol)
end
