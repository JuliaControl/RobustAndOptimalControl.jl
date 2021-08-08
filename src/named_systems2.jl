import ControlSystems as CS
import ControlSystems: nstates, blockdiag

"""
See `named_ss` for a convenient constructor.
"""
struct NamedStateSpace{T,S} <: AbstractStateSpace{T} where S <: AbstractStateSpace{T}
    sys::S
    x
    u
    y
end

function Base.promote_rule(::Type{U}, ::Type{NamedStateSpace{T, S}}) where
    {T, U<:AbstractStateSpace{T} , S<:AbstractStateSpace{T}} 
    inner = promote_type(U,S)
    NamedStateSpace{T, inner}
end

function Base.promote_rule(::Type{NamedStateSpace{T, U}}, ::Type{NamedStateSpace{T, S}}) where
    {T, U<:AbstractStateSpace{T} , S<:AbstractStateSpace{T}} 
    inner = promote_type(U,S)
    NamedStateSpace{T, inner}
end


function Base.convert(::Type{NamedStateSpace{T, S}}, s::U) where {T, S <: AbstractStateSpace, U <: AbstractStateSpace}
    s2 = S === U ? s : Base.convert(S, s)
    named_ss(s2, x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

function Base.convert(::Type{NamedStateSpace{T, S}}, s::NamedStateSpace{T, U}) where {T, S <: AbstractStateSpace, U <: AbstractStateSpace}
    sys = Base.convert(S, s.sys)
    NamedStateSpace{T,typeof(sys)}(sys, s.x, s.u, s.y)
end

# Base.convert(::Type{RobustAndOptimalControl.NamedStateSpace{T, S}}, s::S) where {T, S<:RobustAndOptimalControl.NamedStateSpace} = s



const NamedIndex = Union{Symbol, Vector{Symbol}, Colon}

function Base.getproperty(G::NamedStateSpace, s::Symbol)
    s ∈ fieldnames(NamedStateSpace) && (return getfield(G,s))
    return getproperty(G.sys, s)
end

ControlSystems.numeric_type(G::NamedStateSpace) = ControlSystems.numeric_type(G.sys)

"""
    expand_symbol(s::Symbol, n::Int)

Takes a symbol and an integer and returns a vector of symbols with increasing numbers appended to the end. E.g.,
(:x, 3) -> [:x1, :x2, :x3]

Useful to create signal names for named systems.
"""
expand_symbol(s::Symbol, n::Int) = n == 1 ? [s] : [Symbol(string(s)*string(i)) for i in 1:n]
expand_symbol(v, n::Int) = v

"""
    named_ss(sys::AbstractStateSpace{T}; x, u, y)

Create a `NamedStateSpace` system. This kind of system uses names rather than integer indices to refer to states, inputs and outputs

# Arguments:
- `sys`: A system to add names to.
- `x`: A list of symbols with names of the states.
- `u`: A list of symbols with names of the inputs.
- `y`: A list of symbols with names of the outputs.

Default names of signals if none are provided are `x,u,y`.

# Example
```julia
G1 = ss(1,1,1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
s2 = named_ss(G2, x = :z, u = :u2, y=:y2)

s1[:y1, :u1] # Index using symbols

fb = feedback(s1, s2, r = :r) # 
````
"""
function named_ss(sys::AbstractStateSpace{T};
    x = [Symbol("x$i") for i in 1:sys.nx],
    u = [Symbol("u$i") for i in 1:sys.nu],
    y = [Symbol("y$i") for i in 1:sys.ny],
    ) where T
    x = expand_symbol(x, sys.nx)
    u = expand_symbol(u, sys.nu)
    y = expand_symbol(y, sys.ny)
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
    ii = i isa Colon ? i : names2indices(i, sys.y) 
    jj = j isa Colon ? j : names2indices(j, sys.u) 

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

function Base.:*(s1::Number, s2::NamedStateSpace{T, S}) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        s1*s2.sys,
        s2.x,
        s2.u,
        [Symbol(string(y)*"_scaled") for y in s2.y]
    )
end

function Base.:/(s::NamedStateSpace{T, S}, n::Number) where {T <: CS.TimeEvolution, S}
    (1/n)*s
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
    zi = names2indices(z, P.y)
    yi = names2indices(y, P.y)
    wi = names2indices(w, P.u)
    ui = names2indices(u, P.u)
    ss(P.A, P.B[:, wi], P.B[:, ui], P.C[zi, :], P.C[yi, :], 
        P.D[zi, wi], P.D[zi, ui], P.D[yi, wi], P.D[yi, ui], P.timeevol)
end

function named_ss(sys::ExtendedStateSpace{T};
    x = [Symbol("x$i") for i in 1:sys.nx],
    u = [Symbol("u$i") for i in 1:sys.nu],
    y = [Symbol("y$i") for i in 1:sys.ny],
    w = [Symbol("w$i") for i in 1:sys.nw],
    z = [Symbol("z$i") for i in 1:sys.nz],
    ) where T
    x = expand_symbol(x, sys.nx)
    u = expand_symbol(u, sys.nu)
    y = expand_symbol(y, sys.ny)
    w = expand_symbol(w, sys.nw)
    z = expand_symbol(z, sys.nz)
    length(x) == sys.nx ||
        throw(ArgumentError("Length of state names must match sys.nx ($(sys.nx))"))
    length(u) == sys.nu ||
        throw(ArgumentError("Length of input names must match sys.nu ($(sys.nu))"))
    length(y) == sys.ny ||
        throw(ArgumentError("Length of output names must match sys.ny ($(sys.ny))"))
    length(w) == sys.nw ||
        throw(ArgumentError("Length of disturbance names must match sys.nw ($(sys.nw))"))
    length(z) == sys.nz ||
        throw(ArgumentError("Length of performance names must match sys.nz ($(sys.nz))"))

    sys2 = ss(sys)
    NamedStateSpace{T, typeof(sys2)}(sys2, x, [w; u], [z; y])
end
