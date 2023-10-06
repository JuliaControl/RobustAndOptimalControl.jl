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

import ControlSystemsBase as CS
import ControlSystemsBase: nstates, blockdiag


const NameType = Union{Vector{Symbol}, NTuple{N, Symbol} where N}

"""
See `named_ss` for a convenient constructor.
"""
struct NamedStateSpace{T,S} <: AbstractStateSpace{T} where S <: AbstractStateSpace{T}
    sys::S
    x::Vector{Symbol}
    u::Vector{Symbol}
    y::Vector{Symbol}
    name::String
end

NamedStateSpace(A,B,C,D,x::AbstractVector,u,y,name::String="") = NamedStateSpace{Continuous, StateSpace{Continuous, eltype(A)}}(ss(A,B,C,D), x, u, y, name)
NamedStateSpace(A,B,C,D,Ts::Number,x,u,y,name::String="") = NamedStateSpace{Discrete{typeof(Ts)}, StateSpace{Discrete{typeof(Ts)}, eltype(A)}}(ss(A,B,C,D,Ts), x, u, y,name)

function Base.promote_rule(::Type{U}, ::Type{NamedStateSpace{T, S}}) where
    {T, U<:AbstractStateSpace{T} , S<:AbstractStateSpace{T}} 
    inner = promote_type(U,S)
    NamedStateSpace{T, inner}
end

function Base.promote_rule(::Type{NamedStateSpace{T1, U}}, ::Type{NamedStateSpace{T2, S}}) where
    {T1, T2, U<:AbstractStateSpace{T1} , S<:AbstractStateSpace{T2}} 
    inner = promote_type(U,S)
    TE = promote_type(T1,T2)
    NamedStateSpace{TE, inner}
end

function Base.promote_rule(::Type{NamedStateSpace{TE, StateSpace{TE, T1}}}, ::Type{MT}) where {TE, T1, MT<:AbstractMatrix}
    NamedStateSpace{TE, StateSpace{TE, promote_type(T1,eltype(MT))}}
end



function Base.convert(::Type{NamedStateSpace{T, S}}, s::U) where {T, S <: AbstractStateSpace, U <: AbstractStateSpace}
    s2 = S === U ? s : Base.convert(S, s)
    named_ss(s2, x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

function Base.convert(::Type{NamedStateSpace{T, S}}, M::AbstractMatrix) where {T, S <: AbstractStateSpace}
    te = T <: Discrete ? Discrete(ControlSystemsBase.UNDEF_SAMPLEPETIME) : Continuous()
    NT = numeric_type(S)
    named_ss(ss(NT.(M), te), x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

function Base.convert(::Type{NamedStateSpace{T, S}}, s::NamedStateSpace{T, U}) where {T, S <: AbstractStateSpace, U <: AbstractStateSpace}
    sys = Base.convert(S, s.sys)
    NamedStateSpace{T,typeof(sys)}(sys, s.x, s.u, s.y, s.name)
end

# function Base.convert(::Type{TransferFunction{TE, S}}, s::U) where {TE, S, U <: NamedStateSpace{TE}}
#     convert(TransferFunction{TE, S}, s.sys)
# end

# Base.convert(::Type{RobustAndOptimalControl.NamedStateSpace{T, S}}, s::S) where {T, S<:RobustAndOptimalControl.NamedStateSpace} = s

ControlSystemsBase.system_name(P::NamedStateSpace, i=(:)) = P.name
ControlSystemsBase.input_names(P::NamedStateSpace, i=(:)) = string.(getindex(P.u, i))
ControlSystemsBase.output_names(P::NamedStateSpace, i=(:)) = string.(getindex(P.y, i))
ControlSystemsBase.state_names(P::NamedStateSpace, i=(:)) = string.(getindex(P.x, i))

const NamedIndex = Union{Symbol, Vector{Symbol}, Colon}

function Base.getproperty(G::NamedStateSpace, s::Symbol)
    s ∈ fieldnames(NamedStateSpace) && (return getfield(G,s))
    return getproperty(G.sys, s)
end

Base.propertynames(sys::NamedStateSpace) = (propertynames(sys.sys)..., :x, :u, :y, :name)

ControlSystemsBase.numeric_type(G::NamedStateSpace) = ControlSystemsBase.numeric_type(G.sys)

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

Create a `NamedStateSpace` system. This kind of system uses names rather than integer indices to refer to states, inputs and outputs.
- If a single name is provided but a vector of names is expected, this name will be used as prefix followed by a numerical index.
- If no name is provided, default names (`x,y,u`) will be used.

# Arguments:
- `sys`: A system to add names to.
- `x`: A list of symbols with names of the states.
- `u`: A list of symbols with names of the inputs.
- `y`: A list of symbols with names of the outputs.


# Example
```julia
G1 = ss(1,1,1,0)
G2 = ss(1,1,1,0)
s1 = named_ss(G1, x = :x, u = :u1, y=:y1)
s2 = named_ss(G2, x = :z, u = :u2, y=:y2)

s1[:y1, :u1] # Index using symbols. Uses prefix matching if no exact match is found.

fb = feedback(s1, s2, r = :r) # 
```
"""
function named_ss(sys::AbstractStateSpace{T};
    x = :x,
    u = :u, # sinze is used instead of sys.nu for ExtendedStateSpace
    y = :y,
    name::String = "",
    ) where T
    x = expand_symbol(x, sys.nx)
    u = expand_symbol(u, size(sys,2)) # size is used instead of sys.nu for ExtendedStateSpace
    y = expand_symbol(y, size(sys,1))
    length(x) == sys.nx ||
        throw(ArgumentError("Length of state names must match sys.nx ($(sys.nx)), got length $(length(x))"))
    length(u) == size(sys,2) ||
        throw(ArgumentError("Length of input names must match size(sys,2) ($(size(sys,2))), got length $(length(u))"))
    length(y) == size(sys,1) ||
        throw(ArgumentError("Length of output names must match size(sys,1) ($(size(sys,1))), got length $(length(y))"))

    @check_unique x "x"
    @check_unique u "u"
    @check_unique y "y"

    NamedStateSpace{T, typeof(sys)}(sys, x, u, y, name)
end

"""
    named_ss(sys::AbstractStateSpace, name; x, y, u)

If a single name of the system is provided, the outputs, inputs and states will be automatically named
`y,u,x` with `name` as prefix.
"""
function named_ss(sys::AbstractStateSpace, name;
    x = Symbol(string(name)*"x"),
    y = Symbol(string(name)*"y"),
    u = Symbol(string(name)*"u"),
    )
    named_ss(sys; x, y, u, name=string(name))
end

named_ss(G::LTISystem, args...; kwargs...) = named_ss(ss(G), args...; kwargs...)

ControlSystemsBase.ss(sys::NamedStateSpace) = ss(sys.sys)

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
        sys.name,
    ) # |> sminreal
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
        sys.name,
    ) # |> sminreal
end

function Base.show(io::IO, G::NamedStateSpace)
    name = system_name(G)
    if !isempty(name)
        print(io, name, ": ")
    end
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
        s1.name,
    )
end

function Base.:+(s1::NamedStateSpace{T,S}, s2::NamedStateSpace{T,S}) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        s1.sys+s2.sys,
        [s1.x; s2.x],
        s1.u,
        s1.y,
        "",
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
        "",
    )
end

function Base.:*(s1::Number, s2::NamedStateSpace{T, S}) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        s1*s2.sys,
        s2.x,
        s2.u,
        [Symbol(string(y)*"_scaled") for y in s2.y],
        isempty(s2.name) ? "" : s2.name*"_scaled",
    )
end

function Base.:*(s1::NamedStateSpace{T, S}, s2::Number) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        s1.sys*s2,
        s1.x,
        [Symbol(string(u)*"_scaled") for u in s1.u],
        s1.y,
        isempty(s1.name) ? "" : s1.name*"_scaled",
    )
end

function Base.:/(s::NamedStateSpace{T, S}, n::Number) where {T <: CS.TimeEvolution, S}
    s*(1/n)
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
        "",
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
        "",
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
    sminreal(named_ss(s2; s.x, s.u, y=names, name=s.name))
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
function ControlSystemsBase.feedback(s1::NamedStateSpace{T}, s2::NamedStateSpace{T}; 
    u1=:, w1=:,z1=:,y1=:,u2=:,y2=:,w2=[],z2=[], kwargs...) where {T <: CS.TimeEvolution}
    W1 = names2indices(w1, s1.u)
    U1 = names2indices(u1, s1.u)
    Z1 = names2indices(z1, s1.y)
    Y1 = names2indices(y1, s1.y)

    W2 = names2indices(w2, s2.u)
    U2 = names2indices(u2, s2.u)
    Z2 = names2indices(z2, s2.y)
    Y2 = names2indices(y2, s2.y)

    if s1.nx == 0 && s1.D == I
        # Special case, this happens when calling feedback(I, s2), which calls promote. The promotion creates a named statespace for I with gensym names. This branch gives correct names to the identity such that the sensitivity function will have correct names.
        s1 = named_ss(s1.sys; u = s2.y, y = s2.u, s1.name)
    end

    # TODO: add feedthrough if user requests to have some inputs as outputs

    sys = feedback(s1.sys, s2.sys; W1, W2, U1, U2, Z1, Z2, Y1, Y2, kwargs...)
    fbname = gensym(:feedback)
    x1  = [s1.x; s2.x]
    @check_unique x1
    x1 = [Symbol(string(x1)*string(fbname)) for x1 in x1] # add unique name postfix
    @assert sys.nu == length(W1) + length(W2)
    @assert sys.ny == length(Z1) + length(Z2)
    @assert sys.nx == length(x1)
    nsys = NamedStateSpace{T,typeof(sys)}(sys, x1, s1.u[[W1; W2]], s1.y[[Z1; Z2]], "")
    sminreal(nsys)
end

ControlSystemsBase.feedback(s1::NamedStateSpace{T}, s2::AbstractStateSpace{T}; kwargs...) where {T <: CS.TimeEvolution} = feedback(s1, named_ss(s2); kwargs...)

function connect(systems; u1::Vector{Symbol}, y1::Vector{Symbol}, w1, z1 = (:), verbose = true, kwargs...)
    full = append(systems...)
    @assert length(y1) == length(u1)
    @check_unique u1 "Connected inputs not unique. If you want to connect several signals to the same input, use a summation node, e.g., named_ss(ss([1  1]), u=[:u1, :u2], y=:usum)"
    @check_unique full.u "system inputs"
    @check_unique full.y "system outputs"

    if verbose
        leftover_inputs = setdiff(full.u, [u1; w1])
        isempty(leftover_inputs) || @warn("The following inputs were unconnected $leftover_inputs, ignore this warning if you rely on prefix matching")
        leftover_outputs = setdiff(full.y, z1 == (:) ? y1 : [y1; z1])
        isempty(leftover_outputs) || @warn("The following outputs were unconnected $leftover_outputs, ignore this warning if you rely on prefix matching")
    end


    z2 = []
    w2 = []

    # Connections
    y2 = (:)
    u2 = (:)

    fb = named_ss(ss(I(length(y1)), full.timeevol))
    G = feedback(full, fb; z1, z2, w1, w2, u1, u2, y1, y2, pos_feedback=true, kwargs...)
end


"""
    connect(systems, connections; w1, z1 = (:), verbose = true, kwargs...)

Create block connections using named inputs and outputs.

Addition and subtraction nodes are achieved by creating a linear combination node, i.e., a system with a `D` matrix only.

# Arguments:
- `systems`: A vector of named systems to be connected
- `connections`: a vector of pairs output => input, where each pair maps an output to an input. Each output must appear as an output in one of `systems`, and similarly each input must appear as an input in one of `systems`. All inputs must have unique names and so must all outputs, but an input may have the same name as an output. In the example below the connection `:uP => :uP` connects the output `:uP` of the `addP` block to `P`'s input `:uP`
- `w1`: external signals to be used as inputs in the constructed system. Use `(:)` to indicate all signals
- `z1`: outputs of the constructed system. Use `(:)` to indicate all signals
- `verbose`: Issue warnings for signals that have no connection

Note: Positive feedback is used, controllers that are intended to be connected with negative feedback must thus be negated.

Example:
The following complicated feedback interconnection

```
                 yF
               ┌────────────────────────────────┐
               │                                │
    ┌───────┐  │  ┌───────┐        ┌───────┐    │    ┌───────┐
uF  │       │  │  │       | yR     │       │ yC │ uP │       │   yP
────►   F   ├──┴──►   R   │────┬───►   C   ├────┴────►   P   ├───┬────►
    │       │     │       │    │   │       │         │       │   │
    └───────┘     └───────┘    │   └───────┘         └───────┘   │
                               │                                 │
                               └─────────────────────────────────┘
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
function connect(systems, pairs::AbstractVector{<:Pair}; kwargs...)
    connect(systems; u1 = last.(pairs), y1 = first.(pairs), kwargs...)
end


"""
    splitter(u::Symbol, n::Int, timeevol = Continuous())

Return a named system that splits an input signal into `n` signals. This is useful when an external signal entering a block diagram is to be connected to multiple inputs. See the tutorial 
https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/hinf_connection/
for example usage

# Arguments:
- `u`: Named of the signal to split
- `n`: Number of splits
"""
function splitter(u::Symbol, n::Int, timeevol = Continuous())
    named_ss(ss(ones(n), timeevol), u = [u], y = u^n, name="splitter")
end

# function sumblock(ex::Expr; Ts=0, n=1)
#     @assert ex.head == :(=)
#     sumname = ex.args[1]::Symbol
#     rhs = ex.args[2]::Expr
#     op, s1, s2 = rhs.args
#     timeevol = Ts <= 0 ? ControlSystemsBase.Continuous() : ControlSystemsBase.Discrete(Ts)
#     s = op == :(-) ? -1 : 1
#     named_ss(ss([I(n) s*I(n)], timeevol), u=[s1^n; s2^n], y=sumname)
# end

"""
    sumblock(ex::String; Ts = 0, n = 1)

Create a summation node (named statespace system) that sums (or subtracts) vectors of length `n`.

# Arguments:
- `Ts`: Sample time
- `n`: The length of the input and output vectors. Set `n=1` for scalars.

When using `sumblock` to form block diagrams, note how the system returned from `sumblock` has input names corresponding to the right-hand side of the expression and output names corresponding to the variable on the left-hand side. You will thus typically list connections like `:y => :y` in the connection list to the [`connect`](@ref) function. See the tutorials
- https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/hinf_connection/
- https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/api/#RobustAndOptimalControl.connect-Tuple{Any}
for example usage


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
    timeevol = Ts <= 0 ? ControlSystemsBase.Continuous() : ControlSystemsBase.Discrete(Ts)
    sumname = Symbol(strip(only(match(r"(.+?)\s?=", ex).captures)))^n
    rhs = split(ex, '=', keepempty=false)[2] |> strip
    rhs = replace(rhs, r"([\+\-])" => s" \1 ") # insert whitespace to make sure split works below.
    s = 1
    mats = []
    names = []
    for sym in split(rhs, ' ', keepempty=false)
        if sym == "+"
            s = 1.0
        elseif sym == "-"
            s = -1.0
        else
            push!(mats, s*I(n))
            push!(names, Symbol(sym)^n)
        end
    end
    names = reduce(vcat, names)
    D = reduce(hcat, mats)
    named_ss(ss(D, timeevol), u=names, y=sumname, name="sumblock")
end

function ControlSystemsBase.sminreal(s::NamedStateSpace)
    local sys
    try
        sys = sminreal(s.sys)
    catch
        return s
    end
    _, _, _, inds = CS.struct_ctrb_obsv(s.sys) # we do this one more time to get the inds. This implies repeated calculations, but will allow inner systems of exotic types that have a special method for sminreal to keep their type.
    named_ss(sys; x=s.x[inds], s.u, s.y, s.name)
end

names2indices(::Colon, allnames) = 1:length(allnames) 

function names2indices(names, allnames)
    inds = Union{Nothing, Int}[findfirst(==(n), allnames) for n in names]
    snames = string.(allnames)
    i = 1
    k = 1
    while i <= length(inds)
        if inds[i] === nothing
            # try finding symbols with given prefix
            newi = findall(startswith(string(names[k])), snames)
            newi === nothing || isempty(newi) && error("The indexed NamedSystem has no signal named $(names[i]), available names are $(allnames)")
            deleteat!(inds, i)
            foreach(j->insert!(inds, j+i-1, newi[j]), 1:length(newi))
            i += length(newi)
        else
            i += 1
        end
        k += 1
    end
    inds
end
function names2indices(name::Symbol, allnames)
    i = findfirst(==(name), allnames)
    if i === nothing
        # try finding symbols with given prefix
        i = findall(startswith(string(name)), string.(allnames))
        i === nothing && error("The indexed NamedSystem has no signal named $name, available names are $(allnames)")
        i
    else
        i:i # return a vector rather than scalar for slices of matrices to not drop dim
    end
end

function ExtendedStateSpace(P::NamedStateSpace; z=P.y, y=P.y, w=P.u, u=P.u)
    zi = names2indices(z, P.y)
    yi = names2indices(y, P.y)
    wi = names2indices(w, P.u)
    ui = names2indices(u, P.u)
    ss(P.A, P.B[:, wi], P.B[:, ui], P.C[zi, :], P.C[yi, :], 
        P.D[zi, wi], P.D[zi, ui], P.D[yi, wi], P.D[yi, ui], P.timeevol)
end

"""
    named_ss(sys::ExtendedStateSpace;       kwargs...)
    named_ss(sys::ExtendedStateSpace, name; kwargs...)

Assign names to an ExtendedStateSpace. If no specific names are provided for signals
`z,y,w,u` and states`x`, names will be generated automatically.

# Arguments:
- `name`: Prefix to add to all automatically generated names.
- `x`
- `u`
- `y`
- `w`
- `z`
"""
function named_ss(sys::ExtendedStateSpace{T}, name="";
    x = Symbol(string(name)*"x"),
    u = Symbol(string(name)*"u"),
    y = Symbol(string(name)*"y"),
    w = Symbol(string(name)*"w"),
    z = Symbol(string(name)*"z"),
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
    NamedStateSpace{T, typeof(sys2)}(sys2, x, [w; u], [z; y], name)
end


function partition(P::NamedStateSpace; u=nothing, y=nothing,
    w = nothing,
    z = nothing
)
    if w === nothing
        w = names2indices(setdiff(P.u, u), P.u)
        u = names2indices(u, P.u)
    end
    if z === nothing
        z = names2indices(setdiff(P.y, y), P.y)
        y = names2indices(y, P.y)
    end
    if u === nothing
        u = names2indices(setdiff(P.u, w), P.u)
        w = names2indices(w, P.u)
    end
    if y === nothing
        y = names2indices(setdiff(P.y, z), P.y)
        z = names2indices(z, P.y)
    end
    u = vcat(u)
    y = vcat(y)
    z = vcat(z)
    w = vcat(w)
    ss(P.A, P.B[:, w], P.B[:, u], P.C[z, :], P.C[y, :], 
    P.D[z, w], P.D[z, u], P.D[y, w], P.D[y, u], P.timeevol)
end

function CS.c2d(s::NamedStateSpace{Continuous}, Ts::Real, method::Symbol = :zoh, args...;
    kwargs...)
    named_ss(c2d(s.sys, Ts, method, args...; kwargs...); s.x, s.u, s.y, s.name)
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


function CS.minreal(sys::NamedStateSpace, args...; kwargs...)
    msys = minreal(sys.sys, args...; kwargs...)
    named_ss(msys; sys.u, sys.y, sys.name)
end

for fun in [:baltrunc, :balreal, :balance_statespace]
    @eval function CS.$(fun)(sys::NamedStateSpace, args...; kwargs...)
        msys, rest... = CS.$(fun)(sys.sys, args...; kwargs...)
        named_ss(msys; sys.u, sys.y, sys.name), rest...
    end
end

for fun in [:baltrunc2, :baltrunc_coprime]
    @eval function $(fun)(sys::NamedStateSpace, args...; kwargs...)
        msys, rest... = $(fun)(sys.sys, args...; kwargs...)
        named_ss(msys; sys.u, sys.y, sys.name), rest...
    end
end