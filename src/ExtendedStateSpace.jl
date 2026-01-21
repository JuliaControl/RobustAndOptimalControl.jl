#####################################################################
##                      Data Type Declarations                     ##
#####################################################################
"""
    ExtendedStateSpace{TE, T, S, I} <: AbstractStateSpace{TE}

A type that represents the two-input, two-output system
```
z  ┌─────┐  w
◄──┤     │◄──
   │  P  │
◄──┤     │◄──
y  └─────┘  u
```
where
- `z` denotes controlled outputs (sometimes called performance outputs)
- `y` denotes measured outputs
- `w` denotes external inputs, such as disturbances or references
- `u` denotes control inputs

The call `lft(P, K)` forms the (lower) linear fractional transform
```
z  ┌─────┐  w
◄──┤     │◄──
   │  P  │
┌──┤     │◄─┐
│y └─────┘ u│
│           │
│  ┌─────┐  │
│  │     │  │
└─►│  K  ├──┘
   │     │
   └─────┘
```
i.e., closing the lower loop around `K`.

An `ExtendedStateSpace` can be converted to a standard `StateSpace` by `ss(P)`, this will keep all inputs and outputs, effectively removing the partitioning only.

When [`feedback`](@ref) is called on this type, defaults are automatically set for the feedback indices.
Other functions defined for this type include
- [`system_mapping`](@ref)
- [`performance_mapping`](@ref)
- [`noise_mapping`](@ref)
- [`lft`](@ref)
- [`feedback`](@ref) has special overloads that sets defaults connections for `ExtendedStateSpace`.
and the following design functions expect `ExtendedStateSpace` as inputs
- [`hinfsynthesize`](@ref)
- [`h2synthesize`](@ref)
- [`LQGProblem`](@ref) (also accepts other types)

A video tutorial on how to use this type is available [here](https://youtu.be/huYRrn--AKc).

## Internal Representation
Internally, this type stores an `AbstractStateSpace` along with index vectors that partition the inputs and outputs:
- `sys`: The underlying state-space system
- `w`: Indices for disturbance inputs (corresponds to B1)
- `u`: Indices for control inputs (corresponds to B2)
- `z`: Indices for performance outputs (corresponds to C1)
- `y`: Indices for measured outputs (corresponds to C2)

The type parameter `I` allows index vectors to be `Vector{Int}`, `UnitRange{Int}`, or `Vector{Symbol}` (for NamedStateSpace).
"""
struct ExtendedStateSpace{TE, T, S<:AbstractStateSpace{TE}, I} <: AbstractStateSpace{TE}
    sys::S          # The underlying StateSpace
    w::I            # Disturbance input indices (corresponds to B1)
    u::I            # Control input indices (corresponds to B2)
    z::I            # Performance output indices (corresponds to C1)
    y::I            # Measured output indices (corresponds to C2)
end

# Inner constructor taking individual matrices (for backward compatibility)
function ExtendedStateSpace{TE,T}(
    A,
    B1,
    B2,
    C1,
    C2,
    D11,
    D12,
    D21,
    D22,
    timeevol::TE,
) where {TE,T}
    nx = size(A, 1)
    nw = size(B1, 2)
    nu = size(B2, 2)
    nz = size(C1, 1)
    ny = size(C2, 1)

    size(A, 2) != nx && nx != 0 && error("A must be square")
    size(B1, 1)  == nx ||          error("B1 must have the same row size as A")
    size(B2, 1)  == nx ||          error("B2 must have the same row size as A")
    size(C1, 2)  == nx ||          error("C1 must have the same column size as A")
    size(C2, 2)  == nx ||          error("C2 must have the same column size as A")
    size(D11, 2) == nw ||          error("D11 must have the same column size as B1")
    size(D21, 2) == nw ||          error("D21 must have the same column size as B1")
    size(D12, 2) == nu ||          error("D12 must have the same column size as B2")
    size(D22, 2) == nu ||          error("D22 must have the same column size as B2")
    size(D11, 1) == nz ||          error("D11 must have the same row size as C1")
    size(D12, 1) == nz ||          error("D12 must have the same row size as C1")
    size(D21, 1) == ny ||          error("D21 must have the same row size as C2")
    size(D22, 1) == ny ||          error("D22 must have the same row size as C2")

    # Build combined matrices
    B = [B1 B2]
    C = [C1; C2]
    D = [D11 D12; D21 D22]
    sys = StateSpace{TE, T}(A, B, C, D, timeevol)

    # Compute index vectors
    w_inds = 1:nw
    u_inds = (nw+1):(nw+nu)
    z_inds = 1:nz
    y_inds = (nz+1):(nz+ny)

    ExtendedStateSpace{TE, T, typeof(sys), typeof(w_inds)}(sys, w_inds, u_inds, z_inds, y_inds)
end

# Constructor with all 4 type parameters (for type-preserving operations like negation)
function ExtendedStateSpace{TE, T, S, I}(
    A,
    B1,
    B2,
    C1,
    C2,
    D11,
    D12,
    D21,
    D22,
    timeevol::TE,
) where {TE, T, S, I}
    # Delegate to the simpler constructor
    ExtendedStateSpace{TE, T}(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
end

function ExtendedStateSpace(
    A::AbstractArray{T},
    B1::AbstractArray,
    B2::AbstractArray,
    C1::AbstractArray,
    C2::AbstractArray,
    D11::AbstractArray = zeros(T, size(C1, 1), size(B1, 2)),
    D12::AbstractArray = zeros(T, size(C1, 1), size(B2, 2)),
    D21::AbstractArray = zeros(T, size(C2, 1), size(B1, 2)),
    D22::AbstractArray = zeros(T, size(C2, 1), size(B2, 2)),
    Ts = nothing,
) where T
        # Validate sampling time
    if (Ts isa Real) && Ts <= 0
        error("Ts must be either a positive number or nothing
                (continuous system)")
    end
    if Ts isa Real
        Ts = Discrete(Ts)
    elseif Ts === nothing
        Ts = Continuous()
    end
    return ExtendedStateSpace{typeof(Ts), T}(
        to_matrix(T, A),
        to_matrix(T, B1),
        to_matrix(T, B2),
        to_matrix(T, C1),
        to_matrix(T, C2),
        to_matrix(T, D11),
        to_matrix(T, D12),
        to_matrix(T, D21),
        to_matrix(T, D22),
        Ts,
    )
end

"""
    ExtendedStateSpace(sys::AbstractStateSpace, w, u, z, y)

Create an [`ExtendedStateSpace`](@ref) from an existing state-space system with specified index vectors.

This constructor preserves the type of `sys` (e.g., `NamedStateSpace`).

# Arguments
- `sys`: The underlying state-space system
- `w`: Disturbance input indices (corresponds to B1)
- `u`: Control input indices (corresponds to B2)
- `z`: Performance output indices (corresponds to C1)
- `y`: Measured output indices (corresponds to C2)
"""
function ExtendedStateSpace(sys::S, w::I, u::I, z::I, y::I) where {S<:AbstractStateSpace, I}
    TE = typeof(sys.timeevol)
    T = ControlSystemsBase.numeric_type(sys)
    ExtendedStateSpace{TE, T, S, I}(sys, w, u, z, y)
end

"""
    ss(A, B1, B2, C1, C2, D11, D12, D21, D22 [, Ts])

Create an [`ExtendedStateSpace`](@ref).
"""
function ss(
    A::AbstractArray,
    B1::AbstractArray,
    B2::AbstractArray,
    C1::AbstractArray,
    C2::AbstractArray,
    D11::AbstractArray,
    D12::AbstractArray,
    D21::AbstractArray,
    D22::AbstractArray,
    Ts = nothing,
)
    return ExtendedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22, Ts)
end

function ss(
    A::AbstractArray{T},
    B1::AbstractArray,
    B2::AbstractArray,
    C1::AbstractArray,
    C2::AbstractArray;
    D11::Union{AbstractArray, Real} = 0,
    D12::Union{AbstractArray, Real} = 0,
    D21::Union{AbstractArray, Real} = 0,
    D22::Union{AbstractArray, Real} = 0,
    Ts = nothing,
) where T

    if D11 == 0
        D11 = zeros(T, size(C1, 1), size(B1, 2))
    end
    if D12 == 0
        D12 = zeros(T, size(C1, 1), size(B2, 2))
    end
    if D21 == 0
        D21 = zeros(T, size(C2, 1), size(B1, 2))
    end
    if D22 == 0
        D22 = zeros(T, size(C2, 1), size(B2, 2))
    end
    return ExtendedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22, Ts)
end

function ss(sys::ExtendedStateSpace)
    ss(ssdata(sys)..., sys.timeevol)
end

function Base.promote_rule(::Type{StateSpace{TE, F1}}, ::Type{<:ExtendedStateSpace{TE, F2}}) where {TE, F1, F2}
    ExtendedStateSpace{TE, promote_type(F1, F2)}
end

function Base.convert(::Type{<:ExtendedStateSpace{TE}}, s::StateSpace{TE, F1}) where {TE, F1}
    partition(s, 0, 0)
end

function Base.getproperty(esys::ExtendedStateSpace, s::Symbol)
    # Access to underlying system and index vectors
    if s === :sys || s === :w || s === :u || s === :z || s === :y
        return getfield(esys, s)
    end

    # Get the underlying system for matrix extraction
    sys = getfield(esys, :sys)
    if sys isa NamedStateSpace
        w = names2indices(getfield(esys, :w), sys.u)
        u = names2indices(getfield(esys, :u), sys.u)
        z = names2indices(getfield(esys, :z), sys.y)
        y = names2indices(getfield(esys, :y), sys.y)
    else
        w = getfield(esys, :w)
        u = getfield(esys, :u)
        z = getfield(esys, :z)
        y = getfield(esys, :y)
    end

    # Extract matrices via indexing
    if s === :A
        return sys.A
    elseif s === :B1
        return sys.B[:, w]
    elseif s === :B2
        return sys.B[:, u]
    elseif s === :C1
        return sys.C[z, :]
    elseif s === :C2
        return sys.C[y, :]
    elseif s === :D11
        return sys.D[z, w]
    elseif s === :D12
        return sys.D[z, u]
    elseif s === :D21
        return sys.D[y, w]
    elseif s === :D22
        return sys.D[y, u]
    elseif s === :timeevol
        return sys.timeevol
    elseif s === :Ts
        if isdiscrete(esys)
            return timeevol(esys).Ts
        else
            @warn "Getting time 0.0 for non-discrete systems is deprecated. Check `isdiscrete` before trying to access time."
            return 0.0
        end
    elseif s === :nx
        return nstates(esys)
    elseif s === :nu
        return length(u)
    elseif s === :ny
        return length(y)
    elseif s === :nw
        return length(w)
    elseif s === :nz
        return length(z)
    elseif s === :B
        return [esys.B1 esys.B2]
    elseif s === :C
        return [esys.C1; esys.C2]
    elseif s === :D
        return [esys.D11 esys.D12; esys.D21 esys.D22]
    elseif s === :zinds
        return z
    elseif s === :yinds
        return y
    elseif s === :winds
        return w
    elseif s === :uinds
        return u
    else
        error("type ExtendedStateSpace has no field $s")
    end
end

Base.propertynames(::ExtendedStateSpace) = (:A, :B, :C, :D, :B1, :B2, :C1, :C2, :D11, :D12, :D21, :D22, :Ts, :timeevol, :nx, :ny, :nu, :nw, :nz, :zinds, :yinds, :winds, :uinds, :sys, :w, :u, :z, :y)

ControlSystemsBase.StateSpace(s::ExtendedStateSpace) = ss(ssdata(s)..., s.timeevol)

"""
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(sys)
"""
ssdata_e(sys::ExtendedStateSpace) = sys.A,
sys.B1,
sys.B2,
sys.C1,
sys.C2,
sys.D11,
sys.D12,
sys.D21,
sys.D22

# Funtions for number of intputs, outputs and states
ninputs(sys::ExtendedStateSpace) = size(sys.D, 2)
noutputs(sys::ExtendedStateSpace) = size(sys.D, 1)
nstates(sys::ExtendedStateSpace) = size(sys.A, 1)

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
function Base.:(==)(sys1::ExtendedStateSpace, sys2::ExtendedStateSpace)
    return all(
        getfield(sys1, f) == getfield(sys2, f) for f in fieldnames(ExtendedStateSpace)
    )
end

## Approximate ##
function Base.isapprox(sys1::ExtendedStateSpace, sys2::ExtendedStateSpace)
    return all(
        getfield(sys1, f) ≈ getfield(sys2, f) for f in fieldnames(ExtendedStateSpace)
    )
end

## ADDITION ##
function Base.:+(s1::ExtendedStateSpace, s2::ExtendedStateSpace)
    timeevol = common_timeevol(s1,s2)

    A = blockdiag(s1.A, s2.A)

    B = [[s1.B1; s2.B1] blockdiag(s1.B2, s2.B2)]

    C = [[s1.C1 s2.C1];
    blockdiag(s1.C2, s2.C2)]

    D = [(s1.D11 + s2.D11) s1.D12 s2.D12;
    [s1.D21; s2.D21] blockdiag(s1.D22, s2.D22)]


    P = StateSpace(A, B, C, D, timeevol) # How to handle discrete?
    partition(P, s1.nw, s1.nz)
end

## SUBTRACTION ##
Base.:-(s1::ExtendedStateSpace, s2::ExtendedStateSpace) = +(s1, -s2)

## MULTIPLICATION ##
# NOTE: this series connects the upper parts, i.e., the performance connecmapping
# If you want to series connect the system mapping, use invert_mappings(invert_mappings(sys1)*invert_mappings(sys2))
function Base.:*(s1::ExtendedStateSpace, s2::ExtendedStateSpace)
    timeevol = common_timeevol(s1,s2)

    A = [s1.A                           s1.B1*s2.C1;
    zeros(size(s2.A,1),size(s1.A,2))      s2.A]

    B = [s1.B1*s2.D11                         s1.B2           s1.B1*s2.D12;
    s2.B1              zeros(size(s2.B2,1),size(s1.B2,2))          s2.B2]

    C = [s1.C1                       s1.D11*s2.C1;
    s1.C2                        s1.D21*s2.C1;
    zeros(size(s2.C2,1),size(s1.C2,2))         s2.C2]

    D = [s1.D11*s2.D11           s1.D12        s1.D11*s2.D12;
    s1.D21*s2.D11           s1.D22        s1.D21*s2.D12;
    s2.D21          zeros(size(s2.D22,1),size(s1.D22,2))          s2.D22        ]

    P = StateSpace(A, B, C, D, timeevol)
    partition(P, s2.nw, s1.nz)
end

function Base.:*(s1::ExtendedStateSpace, s2::Number)
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(s1)
    # The reason for only scaling one channel is the use in UncertainSS
    ss(A, s2*B1, B2, C1, C2, s2*D11, D12, s2*D21, D22, s1.timeevol)
    # ExtendedStateSpace(s1.sys*s2, s1.w, s1.u, s1.z, s1.y)
end

function Base.:*(s2::Number, s1::ExtendedStateSpace)
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(s1)
    ss(A, B1, B2, s2*C1, C2, s2*D11, s2*D12, D21, D22, s1.timeevol)
    # ExtendedStateSpace(s2*s1.sys, s1.w, s1.u, s1.z, s1.y)
end

# function invert_mappings(s::ExtendedStateSpace)
#     # Reorder system: inputs [u; w] and outputs [y; z]
#     # This swaps the w↔u and z↔y mappings
#     (; w,u,z,y) = s
#     new_i = [u; w]
#     new_o = [y; z]
#     ExtendedStateSpace(s.sys, s.u, s.w, s.y, s.z)
#     # ExtendedStateSpace(s.sys[new_o, new_i], s.u, s.w, s.y, s.z)
#     # ExtendedStateSpace(s.sys[new_o, new_i], s.w, s.u, s.z, s.y)
# end


function invert_mappings(s::ExtendedStateSpace)
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(s)
    ss(A, B2, B1, C2, C1, D22, D21, D12, D11, s.timeevol)
end

## DIVISION ##
# function Base.:/(n::Number, sys::ExtendedStateSpace)
#     # Ensure s.D is invertible
#     A, B, C, D = ssdata(sys)
#     Dinv = try
#         inv(D)
#     catch
#         error("D isn't invertible")
#     end
#     partition(ss(A - B*Dinv*C, B*Dinv, -n*Dinv*C, n*Dinv, sys.timeevol), sys.nz, sys.nw) # nz, nw reversed (inv)
# end

## NEGATION ##
function Base.:-(sys::ExtendedStateSpace)
    # Negating a StateSpace negates C and D matrices, preserving the internal sys type
    ExtendedStateSpace(-sys.sys, sys.w, sys.u, sys.z, sys.y)
end

#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::ExtendedStateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::ExtendedStateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(sys.D)
Base.size(sys::ExtendedStateSpace, d::Integer) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:ExtendedStateSpace} = S
ControlSystemsBase.numeric_type(sys::ExtendedStateSpace) = eltype(sys.A)

function Base.getindex(sys::ExtendedStateSpace, inds...)
    getindex(sys.sys, inds...)
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

Base.print(io::IO, sys::ExtendedStateSpace) = show(io, sys)

function Base.show(io::IO, sys::ExtendedStateSpace)
    # Compose the name vectors
    println(io, typeof(sys))
    isempty(sys.A) || println(io, "A = \n", _string_mat_with_headers(sys.A))
    isempty(sys.B1) || println(io, "B1 = \n", _string_mat_with_headers(sys.B1))
    isempty(sys.B2) || println(io, "B2 = \n", _string_mat_with_headers(sys.B2))
    isempty(sys.C1) || println(io, "C1 = \n", _string_mat_with_headers(sys.C1))
    isempty(sys.C2) || println(io, "C2 = \n", _string_mat_with_headers(sys.C2))
    isempty(sys.D11) || println(io, "D11 = \n", _string_mat_with_headers(sys.D11))
    isempty(sys.D12) || println(io, "D12 = \n", _string_mat_with_headers(sys.D12))
    isempty(sys.D21) || println(io, "D21 = \n", _string_mat_with_headers(sys.D21))
    isempty(sys.D22) || println(io, "D22 = \n", _string_mat_with_headers(sys.D22))

    # Print sample time
    if isdiscrete(sys) > 0
        println(io, "Sample Time: ", sys.timeevol.Ts, " (seconds)")
    end

    # Print model type
    if nstates(sys) == 0
        print(io, "Static gain")
    elseif iscontinuous(sys)
        print(io, "Continuous-time extended state-space model")
    else
        print(io, "Discrete-time extended state-space model")
    end
end

function ControlSystemsBase.feedback(sys1::ExtendedStateSpace, sys2::AbstractStateSpace;
    W1 = 1:sys1.nw,
    U1 = (1:sys1.nu) .+ sys1.nw,
    Z1 = 1:sys1.nz,
    Y1 = (1:sys1.ny) .+ sys1.nz,
    kwargs...)
    feedback(ss(sys1), sys2; W1, U1, Z1, Y1, kwargs...)
end

function ControlSystemsBase.feedback(sys1::AbstractStateSpace, sys2::ExtendedStateSpace;
    W2 = 1:sys2.nw,
    U2 = (1:sys2.nu) .+ sys2.nw,
    Z2 = 1:sys2.nz,
    Y2 = (1:sys2.ny) .+ sys2.nz,
    kwargs...)
    feedback(sys1, ss(sys2); W2, U2, Z2, Y2, kwargs...)
end

"""
    feedback(sys1::ExtendedStateSpace, sys2::ExtendedStateSpace;
        W1 = 1:sys1.nw,
        U1 = (1:sys1.nu) .+ sys1.nw,
        Z1 = 1:sys1.nz,
        Y1 = (1:sys1.ny) .+ sys1.nz,
        W2 = 1:sys2.nw,
        U2 = (1:sys2.nu) .+ sys2.nw,
        Z2 = 1:sys2.nz,
        Y2 = (1:sys2.ny) .+ sys2.nz,
        kwargs...)

[`ExtendedStateSpace`](@ref) systems use default feedback indices based on the partitioning of the inputs and the outputs.
"""
function ControlSystemsBase.feedback(sys1::ExtendedStateSpace, sys2::ExtendedStateSpace;
    W1 = 1:sys1.nw,
    U1 = (1:sys1.nu) .+ sys1.nw,
    Z1 = 1:sys1.nz,
    Y1 = (1:sys1.ny) .+ sys1.nz,
    W2 = 1:sys2.nw,
    U2 = (1:sys2.nu) .+ sys2.nw,
    Z2 = 1:sys2.nz,
    Y2 = (1:sys2.ny) .+ sys2.nz,
    kwargs...)
    feedback(sys1, ss(sys2); W1, U1, Z1, Y1, W2, U2, Z2, Y2, kwargs...)
end

function ControlSystemsBase.lft(G::ExtendedStateSpace, K, type=:l)

    # if !(G.nu > Δ.ny && G.ny > Δ.nu)
    #     error("Must have G.nu > Δ.ny and G.ny > Δ.nu for lower/upper lft")
    # end

    if K isa ExtendedStateSpace
        error("K isa ExtendedStateSpace not supported")
    end
    if type === :l
        lft(ss(G), K, :l)
    else
        error("Invalid type of lft ($type), specify type=:l")
    end
end

macro sizecompat(a,b)
    n1 = string(a)
    n2 = string(b)
    quote
        s1 = size($(esc(a)), 2)
        s2 = size($(esc(b)), 1)
        if !(s1 == s2)
            throw(ArgumentError(string("size(", $n1, ", 2) and size(", $n2, ", 1) must be compatible, got ($(s1), $(s2))")))
        end
    end
end


function esswrap(f, sys,  args...; kwargs...)
    Gs = ss(sys)
    ret = f(Gs, args...; kwargs...)
    if ret isa AbstractStateSpace
        Gs = ret
        return partition(Gs, sys.nw, sys.nz)
    else
        Gs = first(ret)
        return (partition(Gs, sys.nw, sys.nz), Base.tail(ret)...)
    end
end

ControlSystemsBase.balance_statespace(G::ExtendedStateSpace, perm::Bool=false; kwargs...) = esswrap(ControlSystemsBase.balance_statespace, G, perm; kwargs...)

ControlSystemsBase.similarity_transform(G::ExtendedStateSpace, args...; kwargs...) = esswrap(ControlSystemsBase.similarity_transform, G, args...; kwargs...)

ControlSystemsBase.balreal(G::ExtendedStateSpace, args...; kwargs...) = esswrap(ControlSystemsBase.balreal, G, args...; kwargs...)
modal_form(G::ExtendedStateSpace, args...; kwargs...) = esswrap(modal_form, G, args...; kwargs...)
schur_form(G::ExtendedStateSpace, args...; kwargs...) = esswrap(schur_form, G, args...; kwargs...)
hess_form(G::ExtendedStateSpace, args...; kwargs...) = esswrap(hess_form, G, args...; kwargs...)



# This version is not correct for the intended usage
# function ControlSystemsBase.feedback(s1::ExtendedStateSpace, s2::ExtendedStateSpace;
#     Wperm=:, Zperm=:)

#     timeevol = common_timeevol(s1,s2)
#     s1_B1 = s1.B1
#     s1_B2 = s1.B2
#     s1_C1 = s1.C1
#     s1_C2 = s1.C2
#     s1_D11 = s1.D11
#     s1_D12 = s1.D12
#     s1_D21 = s1.D21
#     s1_D22 = s1.D22


#     # s2_B1 = s2.B2 # These are reversed
#     # s2_B2 = s2.B1
#     # s2_C1 = s2.C2
#     # s2_C2 = s2.C1
#     # s2_D11 = s2.D22
#     # s2_D12 = s2.D21
#     # s2_D21 = s2.D12
#     # s2_D22 = s2.D11

#     s2_B1 = s2.B1 
#     s2_B2 = s2.B2
#     s2_C1 = s2.C1
#     s2_C2 = s2.C2
#     s2_D11 = s2.D11
#     s2_D12 = s2.D12
#     s2_D21 = s2.D21
#     s2_D22 = s2.D22

#     @sizecompat s1_B2 s2_D22
#     @sizecompat s1_B2 s2_C2
#     @sizecompat s2_D22 s1_D21
#     @sizecompat s2_D22 s1_C2
#     @sizecompat s1_D12 s2_C2
#     @sizecompat s2_D12 s1_D22
#     @sizecompat s1_D22 s2_C2

#     if iszero(s1_D22) || iszero(s2_D22)
#         A = [s1.A + s1_B2*s2_D22*s1_C2        s1_B2*s2_C2;
#                 s2_B2*s1_C2            s2.A + s2_B2*s1_D22*s2_C2]

#         B1 = [
#             s1_B1 + s1_B2*s2_D22*s1_D21
#                     s2_B2*s1_D21            
#         ]
#         B2  = [
#             s1_B2*s2_D21
#             s2_B1 + s2_B2*s1_D22*s2_D21
#         ]
#         C1 = [s1_C1+s1_D12*s2_D22*s1_C2        s1_D12*s2_C2]
#         C2 = [s2_D12*s1_C2           s2_C1+s2_D12*s1_D22*s2_C2]
#         D11 = s1_D11 + s1_D12*s2_D22*s1_D21 
#         D12 = s1_D12*s2_D21
#         D21 = s2_D12*s1_D21           
#         D22 = s2_D11 + s2_D12*s1_D22*s2_D21
#     else
#         R1 = lu!(I - s2_D22*s1_D22)
#         issuccess(R1) || 
#             error("Ill-posed feedback interconnection,  I - s2_D22*s1_D22 or I - s2_D22*s1_D22 not invertible")

#         R2 = lu!(I - s1_D22*s2_D22)
#         issuccess(R2) || 
#             error("Ill-posed feedback interconnection,  I - s2_D22*s1_D22 or I - s2_D22*s1_D22 not invertible")

#         A = [s1.A + s1_B2*(R1\s2_D22)*s1_C2        s1_B2*(R1\s2_C2);
#                 s2_B2*(R2\s1_C2)            s2.A + s2_B2*(R2\s1_D22)*s2_C2]

#         B1 = [s1_B1+s1_B2*(R1\s2_D22)*s1_D21;
#                     s2_B2*(R2\s1_D21)]
#         B2 = [s1_B2*(R1\s2_D21); s2_B1 + s2_B2*(R2\s1_D22)*s2_D21]
#         C1 = [s1_C1 + s1_D12*(R1\s2_D22)*s1_C2        s1_D12*(R1\s2_C2)]
#         C2 = [s2_D12*(R2\s1_C2)           s2_C1+s2_D12*(R2\s1_D22)*s2_C2]
#         D11 = [s1_D11 + s1_D12*(R1\s2_D22)*s1_D21]
#         D12 = s1_D12*(R1\s2_D21)
#         D21 = s2_D12*(R2\s1_D21)
#         D22 = s2_D11 + s2_D12*(R2\s1_D22)*s2_D21
#     end

#     return ExtendedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
# end

"""
    se = ExtendedStateSpace(s::AbstractStateSpace; kwargs...)

The conversion from a regular statespace object to an `ExtendedStateSpace` creates the following system by default
```math
\\begin{bmatrix}
    A & B & B \\\\
    C & D & D \\\\
    C & D & D
\\end{bmatrix}
```
i.e., the system and performance mappings are identical, `system_mapping(se) == performance_mapping(se) == s`.
However, all matrices `B1, B2, C1, C2; D11, D12, D21, D22` are overridable by a corresponding keyword argument. In this case, the controlled outputs are the same as measured outputs.

Related: `se = convert(ExtendedStateSpace, s::StateSpace)` produces an `ExtendedStateSpace` with empty `performance_mapping` from w->z such that `ss(se) == s`. `ExtendedStateSpace(sys, B1=I, C1=I)` leads to a system where all state variables are affected by noise, and all are considered performance outputs, this corresponds to the traditional use of the functions `lqr` and `kalman`.
"""
function ExtendedStateSpace(s::AbstractStateSpace;
    A = s.A,
    B1 = s.B,
    B2 = s.B,
    C1 = s.C,
    C2 = s.C,
    D11 = s.D,
    D12 = s.D,
    D21 = s.D,
    D22 = s.D,
    kwargs...
    )
    B1,B2,C1,C2 = _I2mat.((B1,B2,C1,C2), s.nx) # Converts any inputs that are I to appropriate sized matrix
    if size(D11) != (size(C1, 1), size(B1, 2)) && D11 == s.D
        # user provided updated B1 or C1 matrix and this matrix needs to change
        D11 = 0
    end
    if size(D22) != (size(C2, 1), size(B2, 2)) && D22 == s.D
        D22 = 0
    end
    if size(D12) != (size(C1, 1), size(B2, 2)) && D12 == s.D
        # user provided updated B1 or C1 matrix and this matrix needs to change
        D12 = 0
    end
    if size(D21) != (size(C2, 1), size(B1, 2)) && D21 == s.D
        D21 = 0
    end
    ss(copy(A), B1, B2, C1, C2; D11, D12, D21, D22, Ts = s.timeevol, kwargs...)
end
_I2mat(M,nx) = M
_I2mat(i::UniformScaling,nx) = i(nx)

"""
    ExtendedStateSpace(s::ExtendedStateSpace; A, B1, B2, C1, C2, D11, D12, D21, D22, kwargs...)

Create an [`ExtendedStateSpace`](@ref) from an existing one, with the option to override the matrices.
"""
function ExtendedStateSpace(s::ExtendedStateSpace;
    A = s.A,
    B1 = s.B1,
    B2 = s.B2,
    C1 = s.C1,
    C2 = s.C2,
    D11 = s.D11,
    D12 = s.D12,
    D21 = s.D21,
    D22 = s.D22,
    kwargs...
)
    ss(A, B1, B2, C1, C2; D11, D12, D21, D22, Ts = s.timeevol, kwargs...)
end

# ```math
# \\begin{bmatrix}
#     A & I & B \\\\
#     I & 0 & 0 \\\\
#     C & 0 % D
# \\end{bmatrix}
# ```
# i.e., the `system_mapping(se) == s`, while `performance_mapping(s)` 
# ExtendedStateSpace(s::AbstractStateSpace) = ss(s.A, I(s.nx), s.B, I(s.nx), s.C; D22=s.D, Ts = s.timeevol)

"""
    partition(P::AbstractStateSpace; u, y, w=!u, z=!y)

Partition `P` into an [`ExtendedStateSpace`](@ref).
- `u` indicates the indices of the controllable inputs.
- `y` indicates the indices of the measurable outputs.
- `w` indicates the indices of the disturbance inputs (uncontrollable), by default `w` is the complement of `u`.
- `z` indicates the indices of the performance outputs (not neccesarily measurable), by default `z` is the complement of `y`.
"""
function partition(P::AbstractStateSpace; u=nothing, y=nothing,
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
        w = setdiff(1:P.nu, u)
    end
    if z === nothing
        z = setdiff(1:P.ny, y)
    end
    if u === nothing
        u = setdiff(1:P.nu, w)
    end
    if y === nothing
        y = setdiff(1:P.ny, z)
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
    P.timeevol)
end

"""
    partition(P::AbstractStateSpace, nw::Int, nz::Int)
"""
partition(P::AbstractStateSpace, nu1::Int, ny1::Int) = partition(P, w=1:nu1, z=1:ny1)


"""
    system_mapping(P::ExtendedStateSpace)

Return the system from u -> y
See also [`performance_mapping`](@ref), [`system_mapping`](@ref), [`noise_mapping`](@ref)
"""
function system_mapping(P::ExtendedStateSpace, sminreal=sminreal)
    sminreal(P.sys[P.y, P.u])
end

"""
    performance_mapping(P::ExtendedStateSpace)

Return the system from w -> z
See also [`performance_mapping`](@ref), [`system_mapping`](@ref), [`noise_mapping`](@ref)
"""
function performance_mapping(P::ExtendedStateSpace, sminreal=sminreal)
    sminreal(P.sys[P.z, P.w])
end

"""
    noise_mapping(P::ExtendedStateSpace)

Return the system from w -> y
See also [`performance_mapping`](@ref), [`system_mapping`](@ref), [`noise_mapping`](@ref)
"""
function noise_mapping(P::ExtendedStateSpace, sminreal=sminreal)
    sminreal(P.sys[P.y, P.w])
end
