function check_unique(vals, s, msg=""; throw = true)
    u = unique(vals)
    if length(u) != length(vals)
        rep = Dict{Symbol, Int}()
        for v in vals
            n = get(rep, v, 0) + 1
            rep[v] = n
        end
        rep = filter(((_,n),)-> n > 1, pairs(rep))
        repk = keys(rep)
        if throw
            @__MODULE__().throw(ArgumentError(s*" names not unique. Repeated names: "*string(repk)*" "*msg))
        else
            return false
        end
    end
    return true
end

function check_all_unique(s1, s2; throw=true)
    valx = check_unique([getproperty(s1, :x); getproperty(s2, :x)], "x"; throw)
    check_unique([getproperty(s1, :u); getproperty(s2, :u)], "u"; throw=true)
    check_unique([getproperty(s1, :y); getproperty(s2, :y)], "y"; throw=true)
    valx
end

function generate_unique_x_names(systems...)
    x_names = reduce(vcat, s.x for s in systems)
    uniq = check_unique(x_names, "x", throw = false)
    if uniq
        return x_names
    else
        # For each system, if it has a name, append the system name to the x names. If neither system has a name, append gensym names to both systems' x names.
        systemnames = [s.name for s in systems]
        x_names = if count(isempty(s.name) for s in systems) > 1 || length(unique(systemnames)) < length(systemnames)
            # More than one system has empty name or more than one system has the same name
            # We can handle one of the names being empty, which is a common case when the plant has a name but a controller/filter is auto promoted to a named system.
            [gensym(string(x)) for x in x_names]
        else
            # reduce(vcat, [[Symbol(string(s.name)*string(x)) for x in s.x] for s in systems])
            [Symbol(string(s.name)*string(x)) for s in systems for x in s.x]
        end
        @assert allunique(x_names)
        x_names
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
    extra::Dict{Symbol, Any}
end

NamedStateSpace(A,B,C,D,x::AbstractVector,u,y,name::String="",extra=Dict{Symbol,Any}()) = NamedStateSpace{Continuous, StateSpace{Continuous, eltype(A)}}(ss(A,B,C,D), x, u, y, name, extra)
NamedStateSpace(A,B,C,D,Ts::Number,x,u,y,name::String="",extra=Dict{Symbol,Any}()) = NamedStateSpace{Discrete{typeof(Ts)}, StateSpace{Discrete{typeof(Ts)}, eltype(A)}}(ss(A,B,C,D,Ts), x, u, y,name,extra)

# This method is used by the basetype(ST)(A, B, C, D, timeevol) construct
NamedStateSpace(A,B,C,D,te::ControlSystemsBase.TimeEvolution, args...; kwargs...) = named_ss(ss(A,B,C,D,te), args...; kwargs...)

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

function Base.promote_rule(::Type{U}, ::Type{NamedStateSpace{T, S}}) where
    {T, TF, U<:TransferFunction{<:Any, TF} , S<:AbstractStateSpace{T}} 
    inner = promote_type(U,S)
    NamedStateSpace{T, inner}
end

function Base.promote_rule(::Type{NamedStateSpace{TE, StateSpace{TE, T1}}}, ::Type{N}) where {TE, T1, N<:Number}
    NamedStateSpace{TE, StateSpace{TE, promote_type(T1,N)}}
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
    NamedStateSpace{T,typeof(sys)}(sys, s.x, s.u, s.y, s.name, s.extra)
end

function Base.convert(::Type{NamedStateSpace{T, S}}, s::U) where {T, S <: AbstractStateSpace, U <: TransferFunction}
    s2 = Base.convert(S, s)
    named_ss(s2, x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

function Base.convert(::Type{NamedStateSpace{T, S}}, N::Number) where {T, S <: AbstractStateSpace}
    te = T <: Discrete ? Discrete(ControlSystemsBase.UNDEF_SAMPLEPETIME) : Continuous()
    NT = numeric_type(S)
    named_ss(tf(NT(N), te), x = gensym("x"), u = gensym("u"), y = gensym("y"))
end

# function Base.convert(::Type{TransferFunction{TE, S}}, s::U) where {TE, S, U <: NamedStateSpace{TE}}
#     convert(TransferFunction{TE, S}, s.sys)
# end

# Base.convert(::Type{RobustAndOptimalControl.NamedStateSpace{T, S}}, s::S) where {T, S<:RobustAndOptimalControl.NamedStateSpace} = s

ControlSystemsBase.system_name(P::NamedStateSpace, i=(:)) = P.name
ControlSystemsBase.input_names(P::NamedStateSpace, i=(:)) = string.(getindex(P.u, i))
ControlSystemsBase.output_names(P::NamedStateSpace, i=(:)) = string.(getindex(P.y, i))
ControlSystemsBase.state_names(P::NamedStateSpace, i=(:)) = string.(getindex(P.x, i))
function get_extra(P::NamedStateSpace, key::Symbol, default)
    get(P.extra, key, default)
end
function set_extra!(P::NamedStateSpace, key::Symbol, value)
    P.extra[key] = value
    return P
end


"""
    operating_point(P::NamedStateSpace)

Return a named tuple `(; x, u)` containing the operating point of the system `P`. If no operating point is set, a zero operating point of correct dimension is returned.
"""
operating_point(P::NamedStateSpace) = get_extra(P, :operating_point, (x=zeros(P.nx), u=zeros(P.nu)))
operating_point(P::AbstractStateSpace) = (x=zeros(P.nx), u=zeros(P.nu))
function set_operating_point!(P::NamedStateSpace, xu::NamedTuple{(:x,:u)})
    length(xu.x) == P.nx ||
        throw(ArgumentError("Operating point x must have length $(P.nx), got length $(length(xu.x))"))
    length(xu.u) == P.nu ||
        throw(ArgumentError("Operating point u must have length $(P.nu), got length $(length(xu.u))"))
    set_extra!(P, :operating_point, xu)
end
has_operating_point(P::NamedStateSpace) = haskey(P.extra, :operating_point)

"""
    merge_ops(systems...)

Concatenate the operating points of the systems in `systems` into a single operating point. If any system has an operating point, the resulting system will have an operating point with the concatenated state and input vectors.
"""
function merge_ops(systems...)
    if any(has_operating_point, systems)
        opx = reduce(vcat, operating_point(sys).x for sys in systems)
        opu = reduce(vcat, operating_point(sys).u for sys in systems)
        op = (; x = opx, u = opu)
        extra = Dict(:operating_point => op)
    else
        extra = Dict{Symbol, Any}()
    end
end

"""
    merge_ops_x(systems...)

Concatenate the operating points of the systems in `systems` into a single operating point, but only for the state vector `x`. The input vector `u` is taken from the first system's operating point, and all other systems are verified to have the same `u` operating point. If any system has an operating point, the resulting system will have an operating point with the concatenated state vector and the input vector from the first system.
"""
function merge_ops_x(systems...; u_from = :same)
    if any(has_operating_point, systems)
        if u_from === :same
            allequal(operating_point(sys).u for sys in systems) || 
                throw(ArgumentError("All systems must have the same input operating point u to be concatenated."))
            opu = operating_point(systems[1]).u
        elseif u_from isa Integer
            opu = operating_point(systems[u_from]).u
        end
        opx = reduce(vcat, operating_point(sys).x for sys in systems)
        op = (; x = opx, u = opu)
        extra = Dict(:operating_point => op)
    else
        extra = Dict{Symbol, Any}()
    end
end


"""
    infer_operating_point(P1, P2, method = :obsv)

Return the operating point of `P1` inferred from the operating point of `P2` and a similarity transform between the two systems. The method for finding the similarity transform can be specified, default is `:obsv`.

```
s1 = named_ss(ssrand(1,1,2), "P", x = :x, u = :u1, y=:y1)
opx = randn(s1.nx)
opu = randn(s1.nu)
RobustAndOptimalControl.set_operating_point!(s1, (x = opx, u = opu))

s2, _ = balance_statespace(s1) # s1 and s2 are similar systems
op2 = RobustAndOptimalControl.infer_operating_point(s2, s1)
@test RobustAndOptimalControl.operating_point(s2) == op2 
```
"""
function infer_operating_point(P1, P2, method = :obsv)
    T = ControlSystemsBase.find_similarity_transform(P1, P2, method)
    op2 = operating_point(P2)
    op1 = (x = T*op2.x, u = op2.u)
end

function infer_operating_point!(P1, P2, method = :obsv)
    op = infer_operating_point(P1, P2, method)
    set_operating_point!(P1, op)
end

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
    unique = true,
    extra = Dict{Symbol, Any}(),
    ) where T
    if sys isa NamedStateSpace
        error("Cannot wrap a named statespace in a named statespace")
    end
    x = expand_symbol(x, sys.nx)
    u = expand_symbol(u, size(sys,2)) # size is used instead of sys.nu for ExtendedStateSpace
    y = expand_symbol(y, size(sys,1))
    length(x) == sys.nx ||
        throw(ArgumentError("Length of state names must match sys.nx ($(sys.nx)), got length $(length(x))"))
    length(u) == size(sys,2) ||
        throw(ArgumentError("Length of input names must match size(sys,2) ($(size(sys,2))), got length $(length(u))"))
    length(y) == size(sys,1) ||
        throw(ArgumentError("Length of output names must match size(sys,1) ($(size(sys,1))), got length $(length(y))"))

    check_unique(x, "x", "Cannot create a NamedStateSpace system with more than one variable with the same name")
    if unique
        check_unique(y, "y", "To allow connecting a single output signal to several outputs with the same name, pass `unique = false`.")
        check_unique(u, "u", "To allow connecting a single input signal to several inputs with the same name, pass `unique = false`.")
    end

    NamedStateSpace{T, typeof(sys)}(sys, x, u, y, name, extra)
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
    extra = Dict{Symbol, Any}(),
    )
    named_ss(sys; x, y, u, name=string(name), extra)
end

named_ss(G::LTISystem, args...; kwargs...) = named_ss(ss(G), args...; kwargs...)

ControlSystemsBase.ss(sys::NamedStateSpace) = ss(sys.sys)

function ControlSystemsBase.StaticStateSpace(sys::NamedStateSpace) 
    ssys = StaticStateSpace(sys.sys)
    named_ss(ssys, sys.name; sys.x, sys.u, sys.y)
end

iterable(s::Symbol) = [s]
iterable(v) = v

function Base.getindex(sys::NamedStateSpace{T,S}, i::NamedIndex, j::NamedIndex) where {T,S}
    i,j = iterable.((i, j))
    ii = i isa Colon ? i : names2indices(i, sys.y) 
    jj = j isa Colon ? j : names2indices(j, sys.u)

    newsys = NamedStateSpace{T,S}(
        sys.sys[ii, jj],
        sys.x,
        sys.u[jj],
        sys.y[ii],
        sys.name,
        copy(sys.extra),
    ) # |> sminreal
    if has_operating_point(sys)
        op = operating_point(sys)
        op = (; op.x, u=op.u[jj])
        set_operating_point!(newsys, op)
    end
    return newsys
end

function Base.getindex(sys::NamedStateSpace{T,S}, inds...) where {T,S}
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = CS.index2range(inds...)
    newsys = NamedStateSpace{T,S}(
        sys.sys[rows, cols],
        sys.x,
        sys.u[cols],
        sys.y[rows],
        sys.name,
        copy(sys.extra),
    ) # |> sminreal
    if has_operating_point(sys)
        op = operating_point(sys)
        op = (; op.x, u=op.u[cols])
        set_operating_point!(newsys, op)
    end
    return newsys
end

function Base.show(io::IO, G::NamedStateSpace)
    name = system_name(G)
    if !isempty(name)
        print(io, name, ": ")
    end
    print(io, "Named")
    show(io, G.sys)
    length(G.x) < 50 && (print(io, "\nWith state  names: "); println(io, join(G.x, ' ')))
    length(G.u) < 50 && (print(io, "     input  names: "); println(io, join(G.u, ' ')))
    length(G.y) < 50 && (print(io, "     output names: "); println(io, join(G.y, ' ')))
    if has_operating_point(G)
        op = operating_point(G)
        println(io, "Operating point: x = ", op.x, ", u = ", op.u)
    end
end

function Base.:-(s1::NamedStateSpace{T,S}) where {T <: CS.TimeEvolution, S}
    return NamedStateSpace{T,S}(
        -s1.sys,
        s1.x,
        s1.u,
        s1.y,
        s1.name,
        copy(s1.extra),
    )
end

function Base.:+(s1::NamedStateSpace{T,S}, s2::NamedStateSpace{T,S}) where {T <: CS.TimeEvolution, S}
    extra = merge_ops_x(s1, s2)
    return NamedStateSpace{T,S}(
        s1.sys+s2.sys,
        [s1.x; s2.x],
        s1.u,
        s1.y,
        "",
        extra,
    )
end

function Base.:*(s1::NamedStateSpace{T}, s2::NamedStateSpace{T}) where {T <: CS.TimeEvolution}
    x_names = generate_unique_x_names(s1, s2)
    if s1.u != s2.y
        connection_map = join(["$y -> $u" for (u,y) in zip(s1.u, s2.y) if u != y], '\n')
        @warn "Connected signals have different names\n $connection_map" maxlog=2
    end
    sys = s1.sys*s2.sys
    S = typeof(sys)
    extra = merge_ops_x(s1, s2; u_from = 2)
    return NamedStateSpace{T,S}(
        sys,
        x_names,
        s2.u,
        s1.y,
        "",
        extra,
    )
end

function Base.:*(s1::Number, s2::NamedStateSpace{T, S}) where {T <: CS.TimeEvolution, S}
    s3 = s1*s2.sys
    return NamedStateSpace{T, typeof(s3)}(
        s3,
        s2.x,
        s2.u,
        [Symbol(string(y)*"_scaled") for y in s2.y],
        isempty(s2.name) ? "" : s2.name*"_scaled",
        copy(s2.extra),
    )
end

function Base.:*(s1::NamedStateSpace{T, S}, s2::Number) where {T <: CS.TimeEvolution, S}
    s3 = s1.sys*s2
    return NamedStateSpace{T,typeof(s3)}(
        s3,
        s1.x,
        [Symbol(string(u)*"_scaled") for u in s1.u],
        s1.y,
        isempty(s1.name) ? "" : s1.name*"_scaled",
        copy(s1.extra),
    )
end

function Base.:*(s1::AbstractMatrix, s2::NamedStateSpace{T, S}) where {T <: CS.TimeEvolution, S}
    if isdiag(s1)
        return NamedStateSpace{T,S}(
            s1*s2.sys,
            s2.x,
            s2.u,
            [Symbol(string(y)*"_scaled") for y in s2.y],
            isempty(s2.name) ? "" : s2.name*"_scaled",
            copy(s2.extra),
        )
    else
        return *(promote(s1, s2)...)
    end
end

function Base.:*(s1::NamedStateSpace{T, S}, s2::AbstractMatrix) where {T <: CS.TimeEvolution, S}
    if isdiag(s2)
        return NamedStateSpace{T,S}(
            s1.sys*s2,
            s1.x,
            [Symbol(string(u)*"_scaled") for u in s1.u],
            s1.y,
            isempty(s1.name) ? "" : s1.name*"_scaled",
            copy(s1.extra),
        )
    else
        return *(promote(s1, s2)...)
    end
end

function Base.:/(s::NamedStateSpace{T, S}, n::Number) where {T <: CS.TimeEvolution, S}
    s*(1/n)
end

function Base.:/(n::Number, sys::NamedStateSpace)
    isys = n / sys.sys
    if isone(n)
        return named_ss(isys, sys.name*"_inverse", x = sys.x, u = sys.y, y = sys.u)
    else
        return named_ss(isys, sys.name*"_scaled_inverse")
    end
end
##

function Base.hcat(systems::NamedStateSpace{T,S}...) where {T,S}
    x = generate_unique_x_names(systems...)
    u = reduce(vcat, getproperty.(systems, :u))
    check_unique(u, "u")
    extra = merge_ops(systems...)
    return NamedStateSpace{T,S}(
        hcat(getproperty.(systems, :sys)...),
        x,
        u,
        systems[1].y,
        "",
        extra
    )
end

function Base.vcat(systems::NamedStateSpace{T,S}...) where {T,S}
    x = generate_unique_x_names(systems...)
    y = reduce(vcat, getproperty.(systems, :y))
    check_unique(y, "y")
    extra = merge_ops_x(systems...)
    return NamedStateSpace{T,S}(
        vcat(getproperty.(systems, :sys)...),
        x,
        systems[1].u,
        y,
        "",
        extra,
    )
end

"""
    measure(s::NamedStateSpace, names)

Return a system with specified state variables as measurement outputs.

See also [`ControlSystemsBase.add_output`](@ref).
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
    sminreal(named_ss(s2; s.x, s.u, y=names, name=s.name, extra=copy(s.extra)))
end

"""
    merge_nonunique_inputs(sys)

Take a system where one or more input names are appearing more than once, and return a system where the columns of the ``B`` and ``D` matrices corresponding to multiple repeated inputs have been summed together. The resulting system have unique input signal names.

To avoid accidental misuse, a warning is issued if two added colums contain non-zero entries in the same row.

# Example

A system with ``B`` matrix and input names given by
```
B = [1 0 0
     0 2 0
     0 0 3]

u = [:u1, :u2, :u1] # Input names
```
where the input name `:u1` appears more than once, will be reduced to
```
B = [1 0
     0 2
     3 0]
u = [:u1, :u2]
```
"""
function merge_nonunique_inputs(sys)
    i = 0
    inputnames = copy(sys.u)
    while i < length(inputnames)
        i += 1
        inds = findall(n == inputnames[i] for n in inputnames)
        length(inds) == 1 && continue
        # Check that the B-matrix entries are non-overlapping
        Bi = sys.B[:, inds]
        Di = sys.D[:, inds]
        any(>(1), sum(.! iszero.(Bi), dims=2)) && @warn("Input names are not unique and the multiple B-matrix columns associated with the name $(u[i]) have a non-empty intersection of non-zero entries.")
        any(>(1), sum(.! iszero.(Di), dims=2)) && @warn("Input names are not unique and the multiple D-matrix columns associated with the name $(u[i]) have a non-empty intersection of non-zero entries.")
        B = copy(sys.B)
        D = copy(sys.D)
        B[:, inds[1]] = sum(Bi, dims=2)
        D[:, inds[1]] = sum(Di, dims=2)
        kept_inds = setdiff(1:size(B, 2), inds[2:end])
        B = B[:, kept_inds]
        D = D[:, kept_inds]
        deleteat!(inputnames, inds[2:end])
        sys = named_ss(ControlSystemsBase.basetype(sys.sys)(sys.A, B, sys.C, D, sys.timeevol); x=sys.x, u=inputnames, y=sys.y, name=sys.name, unique=false)
    end
    sys
end

"""
    feedback(s1::NamedStateSpace, s2::NamedStateSpace;
            w1 = [],  u1 = (:), z1 = (:), y1 = (:),
            w2 = (:), u2 = (:), z2 = (:), y2 = [], kwargs...)

Feedback between two named systems. The lists of signals to connect are expected to be lists of symbols with signal names, `(:)` indicating all signals, or `[]` indicating no signals.

All signal sets `w1,u1,z1,y1,u2,y2,w2,z2` have the same meaning as for the advanced feedback methods for regular systems.
The added signal set `r` is used to optionally provide a new name for the input of the feedback loop.

To simplify creating complicated feedback interconnections, see `connect`.

If not all inputs and outputs are connected, the returned system may not be a minimal realization.
Use `sminreal` (possibly also `minreal`) to simplify the system if only the input-output map is of interest.
"""
function ControlSystemsBase.feedback(s1::NamedStateSpace{T}, s2::NamedStateSpace{T}; 
    u1=:, w1=:,z1=:,y1=:,u2=:,y2=:,w2=[],z2=[], unique = true, kwargs...) where {T <: CS.TimeEvolution}
    if !unique
        s1 = merge_nonunique_inputs(s1)
    end

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
    x1  = generate_unique_x_names(s1, s2) # [s1.x; s2.x]
    x1 = [Symbol(string(x1)*string(fbname)) for x1 in x1] # add unique name postfix
    @assert sys.nu == length(W1) + length(W2)
    @assert sys.ny == length(Z1) + length(Z2)
    @assert sys.nx == length(x1)
    nsys = NamedStateSpace{T,typeof(sys)}(sys, x1, [s1.u[W1]; s2.u[W2]], [s1.y[Z1]; s2.y[Z2]], "", Dict{Symbol, Any}())
    # sminreal(nsys)
end

ControlSystemsBase.feedback(s1::NamedStateSpace{T}, s2::AbstractStateSpace{T}; kwargs...) where {T <: CS.TimeEvolution} = feedback(s1, named_ss(s2); kwargs...)

function connect(systems; u1::Vector{Symbol}, y1::Vector{Symbol}, external_inputs = nothing, w1 = nothing, external_outputs = (:), z1 = nothing, verbose = true, unique = true, kwargs...)
    full = append(systems...; unique)
    @assert length(y1) == length(u1)

    z1 = something(z1, external_outputs)
    w1 = something(external_inputs, w1)
    w1 === nothing && error("The keyword argument `external_inputs` must be provided")
    if unique
        check_unique(u1, "u1", "Connected inputs not unique. If you want to connect several signals to the same input, use a summation node, e.g., named_ss(ss([1  1]), u=[:u1, :u2], y=:usum)")
        check_unique(full.u, "system inputs", "To allow connecting a single input signal to several inputs with the same name, pass `unique = false`.")
    end

    check_unique(full.y, "system outputs")

    if verbose
        leftover_inputs = setdiff(full.u, [u1; w1])
        isempty(leftover_inputs) || @warn("The following inputs were unconnected $leftover_inputs, ignore this warning if you rely on prefix matching. Turn off this warning by passing `verbose = false`.")
        leftover_outputs = setdiff(full.y, z1 == (:) ? y1 : [y1; z1])
        isempty(leftover_outputs) || @warn("The following outputs were unconnected $leftover_outputs, ignore this warning if you rely on prefix matching. Turn off this warning by passing `verbose = false`.")
    end


    z2 = Symbol[]
    w2 = Symbol[]

    # Connections
    y2 = (:)
    u2 = (:)
    T = numeric_type(full)
    fb = named_ss(ss(one(T)*I(length(y1)), full.timeevol))
    feedback(full, fb; z1, z2, w1, w2, u1, u2, y1, y2, pos_feedback=true, unique, kwargs...)
end


"""
    connect(systems, connections; external_inputs, external_outputs = (:), verbose = true, unique = true, kwargs...)

Create block connections using named inputs and outputs.

Addition and subtraction nodes are achieved by creating a linear combination node, i.e., a system with a `D` matrix only.

# Arguments:
- `systems`: A vector of named systems to be connected
- `connections`: a vector of pairs output => input, where each pair maps an output to an input. Each output must appear as an output in one of `systems`, and similarly each input must appear as an input in one of `systems`. All inputs must have unique names and so must all outputs, but an input may have the same name as an output. In the example below the connection `:uP => :uP` connects the output `:uP` of the `addP` block to `P`'s input `:uP`
- `external_inputs`: external signals to be used as inputs in the constructed system. Use `(:)` to indicate all signals
- `external_outputs`: outputs of the constructed system. Use `(:)` to indicate all signals
- `verbose`: Issue warnings for signals that have no connection
- `unique`: If `true`, all input names must be unique. If `false`, a single external input signal may be connected to multiple input ports with the same name.

Note: Positive feedback is used, controllers that are intended to be connected with negative feedback must thus be negated.

Example:
The following complicated feedback interconnection

```
                 yF
               ┌────────────────────────────────┐
               │                                │
    ┌───────┐  │  ┌───────┐        ┌───────┐    │    ┌───────┐
uF  │       │  │  │       | yR     │       │ yC │ uP │       │   yP
────►   F   ├──┴──►   R   │────+───►   C   ├────+────►   P   ├───┬────►
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
external_inputs = [:uF] # External inputs

G = connect([F, R, C, P, addP, addC], connections; external_inputs)
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
for example usage. An alternative way of connecting an external input to several input ports with the same name is to pass `connect(..., unique=false)`.

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
    op = operating_point(s)
    op = (x = op.x[inds], u = op.u)
    newsys = named_ss(sys; x=s.x[inds], s.u, s.y, s.name, extra=copy(s.extra))
    set_operating_point!(newsys, op)
    return newsys
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
    NamedStateSpace{T, typeof(sys2)}(sys2, x, [w; u], [z; y], name, Dict{Symbol, Any}())
end


function partition(P::NamedStateSpace; u=nothing, y=nothing,
    w = nothing,
    z = nothing,
    B1 = nothing,
    B2 = nothing,
    C1 = nothing,
    C2 = nothing,
    D11 = nothing,
    D12 = nothing,
    D21 = nothing,
    D22 = nothing,
)
    if w === nothing
        inds = names2indices(u, P.u)
        w = setdiff(1:P.nu, inds)
        u = inds
    end
    if z === nothing
        inds = names2indices(y, P.y)
        z = setdiff(1:P.ny, inds)
        y = inds
    end
    if u === nothing
        inds = names2indices(w, P.u)
        u = setdiff(1:P.nu, inds)
        w = inds
    end
    if y === nothing
        inds = names2indices(z, P.y)
        y = setdiff(1:P.ny, inds)
        z = inds
    end
    u = vcat(u)
    y = vcat(y)
    z = vcat(z)
    w = vcat(w)
    ss(P.A,
        B1 === nothing ? P.B[:, w] : B1,
        B2 === nothing ? P.B[:, u] : B2,
        C1 === nothing ? P.C[z, :] : C1,
        C2 === nothing ? P.C[y, :] : C2 ,
        D11 === nothing ? P.D[z, w] : D11,
        D12 === nothing ? P.D[z, u] : D12,
        D21 === nothing ? P.D[y, w] : D21,
        D22 === nothing ? P.D[y, u] : D22,
        P.timeevol
    )
end

function CS.c2d(s::NamedStateSpace{Continuous}, Ts::Real, method::Symbol = :zoh, args...;
    kwargs...)
    named_ss(c2d(s.sys, Ts, method, args...; kwargs...); s.x, s.u, s.y, s.name, extra=copy(s.extra))
end


function CS.append(systems::NamedStateSpace...; kwargs...)
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
    
    extra = merge_ops(systems...)
    newsys = named_ss(systype(A, B, C, D, timeevol); x, y, u, extra, kwargs...)
    return newsys
end


"""
    add_output(sys::NamedStateSpace, C2::AbstractArray, D2 = 0; y)

Add outputs to `sys` corresponding to the output matrix `C2` and the feedthrough matrix `D2` to the system `sys`.

# Arguments:
- `y`: The names used for the new outputs. If not provided, the names will be generated automatically.

See also [`measure`](@ref) for a simpler way to output state variables.
"""
function CS.add_output(sys::NamedStateSpace, C2::AbstractArray, D2=0; y = [Symbol("y_$i") for i in (1:size(C2, 1)) .+ sys.ny])
    T = promote_type(CS.numeric_type(sys), eltype(C2), eltype(D2))
    A,B,C,D = ssdata(sys)
    D3 = D2 == 0 ? zeros(T, size(C2, 1), sys.nu) : D2
    x = sys.x
    u = sys.u
    y = [sys.y; y]
    named_ss(ss(A, B, [C; C2], [D; D3]), sys.timeevol; x, u, y, extra=copy(sys.extra))
end


function CS.minreal(sys::NamedStateSpace, args...; kwargs...)
    msys = minreal(sys.sys, args...; kwargs...)
    named_ss(msys; sys.u, sys.y, sys.name)
end

for fun in [:baltrunc, :balreal]
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

function CS.balance_statespace(sys::NamedStateSpace, args...; kwargs...)
    msys, T, rest... = CS.balance_statespace(sys.sys, args...; kwargs...)
    newsys = named_ss(msys; sys.u, sys.y, sys.name, extra=copy(sys.extra))
    if has_operating_point(sys)
        op = operating_point(sys)
        op = (x = T*op.x, u = op.u)
        set_operating_point!(newsys, op)
    end
    newsys, T, rest...
end