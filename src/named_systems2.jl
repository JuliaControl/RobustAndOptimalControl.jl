macro check_unique(ex, s = string(ex), msg="")
    quote
        $(esc(__source__))
        vals = $(esc(ex))
        u = unique(vals)
        if length(u) != length(vals)
            rep = Dict{Symbol, Int}()
            for v in vals
                n = get(rep, v, 0) + 1
                rep[v] = n
            end
            rep = filter(((_,n),)-> n > 1, pairs(rep))
            repk = keys(rep)
            throw(ArgumentError($(s)*" names not unique. Repeated names: "*string(repk)*" "*$(esc(msg))))
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

import ControlSystems as CS
import ControlSystems: nstates, blockdiag


const NameType = Union{Vector{Symbol}, NTuple{N, Symbol} where N}

"""
See `named_ss` for a convenient constructor.
"""
struct NamedStateSpace{T,S} <: AbstractStateSpace{T} where S <: AbstractStateSpace{T}
    sys::S
    x::Vector{Symbol}
    u::Vector{Symbol}
    y::Vector{Symbol}
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

The short-hand syntax `s^n` is also available, e.g., `:x^3 == expand_symbol(:x, 3)`.

Useful to create signal names for named systems.
"""
expand_symbol(s::Symbol, n::Int) = n == 1 ? [s] : [Symbol(string(s)*string(i)) for i in 1:n]
expand_symbol(v, n::Int) = v
Base.:^(s::Symbol, n::Int) = expand_symbol(s, n)

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

    @check_unique x "x"
    @check_unique u "u"
    @check_unique y "y"

    NamedStateSpace{T, typeof(sys)}(sys, x, u, y)
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
    ) |> sminreal
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
    ) |> sminreal
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
        @warn "Connected signals have different names\n $connection_map" maxlog=2
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

"""
    measure(s::NamedStateSpace, names)

Return a system with specified states as measurement outputs.
"""
function measure(s::NamedStateSpace, names)
    inds = names2indices(names, s.x)
    A,B = ssdata(s)
    ny = length(inds)
    C = zeros(eltype(A), ny, s.nx)
    for i = 1:ny
        C[i, inds[i]] = 1
    end
    s2 = ss(A,B,C,0, s.timeevol)
    sminreal(named_ss(s2; s.x, s.u, y=names))
end

"""
    feedback(s1::NamedStateSpace, s2::NamedStateSpace;
            w1 = [],  u1 = (:), z1 = (:), y1 = (:),
            w2 = (:), u2 = (:), z2 = (:), y2 = [], kwargs...)

Feedback between two named systems. The lists of signals to connect are expected to be lists of symbols with signal names, `(:)` indicating all signals, or `[]` indicating no signals.

All signal sets `w1,u1,z1,y1,u2,y2,w2,z2` have the same meaning as for the advanced feedback methods for regular systems.
The added signal set `r` is used to optionally provide a new name for the input of the feedback loop.

To simplify creating complicated feedback interconnections, see `connect`.
"""
function ControlSystems.feedback(s1::NamedStateSpace{T}, s2::NamedStateSpace{T}; 
    u1=:, w1=[],z1=:,y1=:,u2=:,y2=:,w2=:,z2=[], kwargs...) where {T <: CS.TimeEvolution}
    W1 = names2indices(w1, s1.u)
    U1 = names2indices(u1, s1.u)
    Z1 = names2indices(z1, s1.y)
    Y1 = names2indices(y1, s1.y)

    W2 = names2indices(w2, s2.u)
    U2 = names2indices(u2, s2.u)
    Z2 = names2indices(z2, s2.y)
    Y2 = names2indices(y2, s2.y)

    # TODO: add feedthrough if user requests to have some inputs as outputs

    sys = feedback(s1.sys, s2.sys; W1, W2, U1, U2, Z1, Z2, Y1, Y2, kwargs...)
    fbname = gensym(:feedback)
    x1  = [s1.x; s2.x]
    @check_unique x1
    x1 = [Symbol(string(x1)*string(fbname)) for x1 in x1] # add unique name postfix
    @assert sys.nu == length(W1) + length(W2)
    @assert sys.ny == length(Z1) + length(Z2)
    @assert sys.nx == length(x1)
    nsys = NamedStateSpace{T,typeof(sys)}(sys, x1, s1.u[[W1; W2]], s1.y[[Z1; Z2]])
    sminreal(nsys)
end

"""
    connect(systems, connections; w1, z1 = (:), verbose = true, kwargs...)

    Create complicated feedback interconnection. 

Addition and subtraction nodes are achieved by creating a linear combination node, i.e., a system with a `D` matrix only.

# Arguments:
- `systems`: A vector of named systems to be connected
- `connections`: a vector of pairs indicating output => input mappings.
    - `u1`: input mappings  (alternative input argument)
    - `y1`: output mappings (alternative input argument)
- `w1`: external signals
- `z1`: outputs (can overlap with `y1`)
- `verbose`: Issue warnings for signals that have no connection

Example:
The following complicated feedback interconnection

```
                 yF
              ┌────────────────────────────────┐
              │                                │
    ┌───────┐ │  ┌───────┐ yR   ┌─────────┐    │    ┌───────┐
uF  │       │ │  │       ├──────►         │ yC │  uP│       │    yP
────►   F   ├─┴──►   R   │      │    C    ├────+────►   P   ├────┬────►
    │       │    │       │   ┌──►         │         │       │    │
    └───────┘    └───────┘   │  └─────────┘         └───────┘    │
                             │                                   │
                             └───────────────────────────────────┘
```
can be created by
```
F = named_ss(ssrand(1, 1, 2, proper=true), x=:xF, u=:uF, y=:yF)
R = named_ss(ssrand(1, 1, 2, proper=true), x=:xR, u=:uR, y=:yR)
C = named_ss(ssrand(1, 1, 2, proper=true), x=:xC, u=:uC, y=:yC)
P = named_ss(ssrand(1, 1, 3, proper=true), x=:xP, u=:uP, y=:yP)

addP = sumblock("uP = yF + yC") # Sum node before P
addC = sumblock("uC = yR - yP") # Sum node before C

connections = [
    :yP => :yP # Output to input
    :uP => :uP
    :yC => :yC
    :yF => :yF
    :yF => :uR
    :uC => :uC
    :yR => :yR
]
w1 = [:uF] # External inputs

G = connect([F, R, C, P, addP, addC], connections; w1)
```

If an external input is to be connected to multiple points, use a `splitter` to split up the signal into a set of unique names which are then used in the connections.
"""
function connect(systems; u1::Vector{Symbol}, y1::Vector{Symbol}, w1::Vector{Symbol}, z1 = (:), verbose = true, kwargs...)
    full = append(systems...)
    @assert length(y1) == length(u1)
    @check_unique u1 u1 "Connected inputs not unique. If you want to connect several signals to the same input, use a summation node, e.g., named_ss(ss([1  1]), u=[:u1, :u2], y=:usum)"
    @check_unique full.u "system inputs"
    @check_unique full.y "system outputs"

    if verbose
        leftover_inputs = setdiff(full.u, [u1; w1])
        isempty(leftover_inputs) || @warn("The following inputs were unconnected $leftover_inputs")
        leftover_outputs = setdiff(full.y, z1 == (:) ? y1 : [y1; z1])
        isempty(leftover_outputs) || @warn("The following outputs were unconnected $leftover_outputs")
    end
    

    z2 = []
    w2 = []

    # Connections
    y2 = (:)
    u2 = (:)

    fb = named_ss(ss(I(length(y1)), full.timeevol))
    G = feedback(full, fb; z1, z2, w1, w2, u1, u2, y1, y2, pos_feedback=true, kwargs...)
end

function connect(systems, pairs::AbstractVector{<:Pair}; kwargs...)
    connect(systems; u1 = last.(pairs), y1 = first.(pairs), kwargs...)
end

function splitter(u::Symbol, n::Int, timeevol = Continuous())
    named_ss(ss(ones(n), timeevol), u = [u], y = u^n)
end

# function sumblock(ex::Expr; Ts=0, n=1)
#     @assert ex.head == :(=)
#     sumname = ex.args[1]::Symbol
#     rhs = ex.args[2]::Expr
#     op, s1, s2 = rhs.args
#     timeevol = Ts <= 0 ? ControlSystems.Continuous() : ControlSystems.Discrete(Ts)
#     s = op == :(-) ? -1 : 1
#     named_ss(ss([I(n) s*I(n)], timeevol), u=[s1^n; s2^n], y=sumname)
# end

"""
    sumblock(ex::String; Ts = 0, n = 1)

Create a summation node that sums (or subtracts) vectors of length `n`.

# Arguments:
- `Ts`: Sample time
- `n`: The length of the input and output vectors. Set `n=1` for scalars.

# Examples:
```
julia> sumblock("uP = vf + yL")
NamedStateSpace{Continuous, Int64}
D = 
 1  1

With state  names: 
     input  names: vf yL
     output names: uP


julia> sumblock("x_diff = xr - xh"; n=3)
NamedStateSpace{Continuous, Int64}
D = 
 1  0  0  -1   0   0
 0  1  0   0  -1   0
 0  0  1   0   0  -1

With state  names: 
     input  names: xr1 xr2 xr3 xh1 xh2 xh3
     output names: x_diff1 x_diff2 x_diff3
     

julia> sumblock("a = b + c - d")
NamedStateSpace{Continuous, Int64}
D = 
 1  1  -1

With state  names: 
     input  names: b c d
     output names: a
```
"""
function sumblock(ex::String; Ts=0, n=1)
    timeevol = Ts <= 0 ? ControlSystems.Continuous() : ControlSystems.Discrete(Ts)
    sumname = Symbol(strip(only(match(r"(.+?)\s?=", ex).captures)))^n
    rhs = split(ex, '=', keepempty=false)[2] |> strip
    rhs = replace(rhs, r"([\+\-])" => s" \1 ") # insert whitespace to make sure split works below.
    s = 1
    mats = []
    names = []
    for sym in split(rhs, ' ', keepempty=false)
        if sym == "+"
            s = 1
        elseif sym == "-"
            s = -1
        else
            push!(mats, s*I(n))
            push!(names, Symbol(sym)^n)
        end
    end
    names = reduce(vcat, names)
    D = reduce(hcat, mats)
    named_ss(ss(D, timeevol), u=names, y=sumname)
end

function ControlSystems.sminreal(s::NamedStateSpace)
    local sys
    try
        sys = sminreal(s.sys)
    catch
        return s
    end
    _, _, _, inds = CS.struct_ctrb_obsv(s.sys) # we do this one more time to get the inds. This implies repeated calculations, but will allow inner systems of exotic types that have a special method for sminreal to keep their type.
    named_ss(sys; x=s.x[inds], s.u, s.y)
end

names2indices(::Colon, allnames) = 1:length(allnames) 

function names2indices(names, allnames)
    inds = [findfirst(==(n), allnames) for n in names]
    for i in eachindex(inds)
        inds[i] === nothing && error("The indexed NamedSystem has no signal named $(names[i]), available names are $(allnames)")
    end
    inds
end
function names2indices(name::Symbol, allnames)
    i = findfirst(==(name), allnames)
    i === nothing && error("The indexed NamedSystem has no signal named $name, available names are $(allnames)")
    i:i # return a vector rather than scalar for slices of matrices to not drop dim
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

function CS.stepplot(s::NamedStateSpace, args...;
    title  = permutedims(["Step Response from $n" for n in s.u]),
    kwargs...)
    stepplot(s.sys, args...; kwargs...)
    CS.Plots.plot!(;
        title,
        ylabel = permutedims(["$n" for n in s.y]),
    )
end

function CS.lsimplot(s::NamedStateSpace, args...;
    title  = permutedims(["From $n" for n in s.u]),
    kwargs...)
    lsimplot(s.sys, args...; kwargs...)
    CS.Plots.plot!(;
        title,
        ylabel = permutedims(["$n" for n in s.y]),
    )
end

function CS.bodeplot(s::NamedStateSpace, args...;
    title  = permutedims(["From $n" for n in s.u]),
    kwargs...)
    bodeplot(s.sys, args...; kwargs...)
    CS.Plots.plot!(;
        title,
        ylabel = permutedims(["$n" for n in s.y]),
    )
end

function CS.c2d(s::NamedStateSpace, args...;
    kwargs...)
    named_ss(c2d(s.sys, args...; kwargs...); s.x, s.u, s.y)
end


function CS.append(systems::NamedStateSpace...)
    systype = promote_type([typeof(s.sys) for s in systems]...)
    timeevol = common_timeevol(systems...)
    A = blockdiag([s.A for s in systems]...)
    B = blockdiag([_remove_empty_cols(s.B) for s in systems]...)
    C = blockdiag([s.C for s in systems]...)
    D = blockdiag([_remove_empty_cols(s.D) for s in systems]...)

    # TODO: If some name below is repeated, perhaps the B matrix can be reduced in width?
    # This could be optional, i.e., unique_names=true
    x = reduce(vcat, getproperty.(systems, :x))
    y = reduce(vcat, getproperty.(systems, :y))
    u = reduce(vcat, getproperty.(systems, :u))

    return named_ss(systype(A, B, C, D, timeevol); x, y, u)
end
