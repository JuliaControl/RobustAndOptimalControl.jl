import ControlSystems as CS
import ControlSystems: nstates, blockdiag

struct NamedStateSpace{T} <: AbstractStateSpace{T}
    sys::AbstractStateSpace{T}
    x_names
    u_names
    y_names
end
const NamedIndex = Union{Symbol, Vector{Symbol}}

function Base.getproperty(G::NamedStateSpace, s::Symbol)
    s ∈ fieldnames(NamedStateSpace) && (return getfield(G,s))
    return getproperty(G.sys, s)
end

ControlSystems.numeric_type(G::NamedStateSpace) = ControlSystems.numeric_type(G.sys)

function named_ss(sys;
    x_names = [Symbol("x$i") for i in 1:sys.nx],
    u_names = [Symbol("u$i") for i in 1:sys.nu],
    y_names = [Symbol("y$i") for i in 1:sys.ny],
    )
    length(x_names) == sys.nx  || throw(ArgumentError("Length of state names must match sys.nx"))
    length(u_names) == sys.nu  || throw(ArgumentError("Length of input names must match sys.nu"))
    length(y_names) == sys.ny || throw(ArgumentError("Length of state names must match sys.ny"))

    NamedStateSpace(sys, x_names, u_names, y_names)
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
        vals = [getproperty($(esc(s1)), :x_names); getproperty($(esc(s2)), :x_names)]
        @check_unique vals "x_names"
        vals = [getproperty($(esc(s1)), :u_names); getproperty($(esc(s2)), :u_names)]
        @check_unique vals "u_names"
        vals = [getproperty($(esc(s1)), :y_names); getproperty($(esc(s2)), :y_names)]
        @check_unique vals "y_names"
    end
end

iterable(s::Symbol) = [s]
iterable(v) = v

function Base.getindex(sys::NamedStateSpace, i::NamedIndex, j::NamedIndex)
    i,j = iterable.((i, j))
    ii = findall(sys.y_names .∈ (i, )) # findall(i .∈ (sys.y_names, ))
    jj = findall(sys.u_names .∈ (j, )) # findall(j .∈ (sys.u_names, ))

    return NamedStateSpace(
        sys.sys[ii, jj],
        sys.x_names,
        sys.u_names[jj],
        sys.y_names[ii],
    )
end

function Base.show(io::IO, G::NamedStateSpace)
    show(io, G.sys)
    print(io, "\nWith state  names: "); println(io, join(G.x_names, ' '))
    print(io, "     input  names: "); println(io, join(G.u_names, ' '))
    print(io, "     output names: "); println(io, join(G.y_names, ' '))
end

function Base.:-(s1::NamedStateSpace{TE}) where {TE <: CS.TimeEvolution}
    return NamedStateSpace(
        -s1.sys,
        s1.x_names,
        s1.u_names,
        s1.y_names,
    )
end

function Base.:+(s1::NamedStateSpace{TE}, s2::NamedStateSpace{TE}) where {TE <: CS.TimeEvolution}
    return NamedStateSpace(
        s1.sys+s2.sys,
        [s1.x_names; s2.x_names],
        s1.u_names,
        s1.y_names,
    )
end

function Base.:*(s1::NamedStateSpace{TE}, s2::NamedStateSpace{TE}) where {TE <: CS.TimeEvolution}
    @check_all_unique s1 s2
    return NamedStateSpace(
        s1.sys*s2.sys,
        [s1.x_names; s2.x_names],
        s2.u_names,
        s1.y_names,
    )
end
##

function Base.hcat(systems::NamedStateSpace{<:TE}...) where TE
    x_names = reduce(vcat, getproperty.(systems, :x_names))
    u_names = reduce(vcat, getproperty.(systems, :u_names))
    @check_unique x_names
    @check_unique u_names
    return NamedStateSpace(
        hcat(getproperty.(systems, :sys)...),
        x_names,
        u_names,
        systems[1].y_names,
    )
end

function Base.vcat(systems::NamedStateSpace{<:TE}...) where TE
    x_names = reduce(vcat, getproperty.(systems, :x_names))
    y_names = reduce(vcat, getproperty.(systems, :y_names))
    @check_unique x_names
    @check_unique y_names
    return NamedStateSpace(
        vcat(getproperty.(systems, :sys)...),
        x_names,
        systems[1].u_names,
        y_names,
    )
end


function ControlSystems.feedback(s1::NamedStateSpace{TE}, s2::NamedStateSpace{TE}; r_names = s1.u_names) where {TE <: CS.TimeEvolution}
    sys = feedback(s1.sys, s2.sys)
    r_names = iterable(r_names)
    x1_names  = [s1.x_names; s2.x_names]
    @check_unique x1_names
    return NamedStateSpace(
        sys,
        x1_names,
        r_names,
        s1.y_names,
    )
end
