#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

struct ExtendedStateSpace{TE,T} <: AbstractStateSpace{TE}
    A::Matrix{T}
    B1::Matrix{T}
    B2::Matrix{T}
    C1::Matrix{T}
    C2::Matrix{T}
    D11::Matrix{T}
    D12::Matrix{T}
    D21::Matrix{T}
    D22::Matrix{T}
    timeevol::TE
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

        if size(A, 2) != nx && nx != 0
            error("A must be square")
        elseif size(B1, 1) != nx
            error("B1 must have the same row size as A")
        elseif size(B2, 1) != nx
            error("B2 must have the same row size as A")
        elseif size(C1, 2) != nx
            error("C1 must have the same column size as A")
        elseif size(C2, 2) != nx
            error("C2 must have the same column size as A")
        elseif nw != size(D11, 2)
            error("D11 must have the same column size as B1")
        elseif nw != size(D21, 2)
            error("D12 must have the same column size as B1")
        elseif nu != size(D12, 2)
            error("D12 must have the same column size as B2")
        elseif nu != size(D22, 2)
            error("D22 must have the same column size as B2")
        elseif nz != size(D11, 1)
            error("D11 must have the same row size as C1")
        elseif nz != size(D12, 1)
            error("D12 must have the same row size as C1")
        elseif ny != size(D21, 1)
            error("D11 must have the same row size as C12")
        elseif ny != size(D22, 1)
            error("D12 must have the same row size as C2")
        end


        new{TE,T}(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
    end
end

function ExtendedStateSpace(
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
    T = Float64
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

# function ExtendedStateSpace(
#     A::AbstractArray,
#     B1::AbstractArray,
#     B2::AbstractArray,
#     C1::AbstractArray,
#     C2::AbstractArray,
#     D11::AbstractArray,
#     D12::AbstractArray,
#     D21::AbstractArray,
#     D22::AbstractArray,
#     Ts::Real,
# )
#     # TODO: change back in 0.7 T = promote_type(eltype(A),eltype(B),eltype(C),eltype(D))
#     TBC = promote_type(
#         promote_type(eltype(B1), eltype(B2)),
#         promote_type(eltype(C1), eltype(C2)),
#     )
#     TD = promote_type(
#         promote_type(eltype(D11), eltype(D12)),
#         promote_type(eltype(D21), eltype(D22)),
#     )
#     T = promote_type(promote_type(TBC, TD), eltype(A))
#     @assert (
#         typeof(to_matrix(T, A)) ==
#         typeof(to_matrix(T, B1)) ==
#         typeof(to_matrix(T, B2)) ==
#         typeof(to_matrix(T, C1)) ==
#         typeof(to_matrix(T, C2)) ==
#         typeof(to_matrix(T, D11)) ==
#         typeof(to_matrix(T, D12)) ==
#         typeof(to_matrix(T, D21))
#     )
#     return ExtendedStateSpace{T,Matrix{T}}(
#         to_matrix(T, A),
#         to_matrix(T, B1),
#         to_matrix(T, B2),
#         to_matrix(T, C1),
#         to_matrix(T, C2),
#         to_matrix(T, D11),
#         to_matrix(T, D12),
#         to_matrix(T, D21),
#         to_matrix(T, D22),
#         Float64(Ts),
#     )
# end

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


function Base.getproperty(sys::ExtendedStateSpace, s::Symbol)
    if s === :Ts
        # if !isdiscrete(sys) # NOTE this line seems to be breaking inference of isdiscrete (is there a test for this?)
        if isdiscrete(sys)
            return timeevol(sys).Ts
        else
            @warn "Getting time 0.0 for non-discrete systems is deprecated. Check `isdiscrete` before trying to access time."
            return 0.0
        end
    elseif s === :nx
        return nstates(sys)
    elseif s === :nu
        return ninputs(sys)
    elseif s === :ny
        return size(sys.C2, 1)
    elseif s === :nw
        return size(sys.B1, 2)
    elseif s === :nz
        return size(sys.C1, 1)
    elseif s === :B
        [sys.B1 sys.B2]
    elseif s === :C
        [sys.C1; sys.C2]
    elseif s === :D
        [sys.D11 sys.D12; sys.D21 sys.D22]
    else
        return getfield(sys, s)
    end
end

ControlSystems.StateSpace(s::ExtendedStateSpace) = ss(ssdata(s)...)

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
# not sure how to best handle this yet

## SUBTRACTION ##
# not sure how to best handle this yet

## NEGATION ##
# not sure how to best handle this yet

## MULTIPLICATION ##
# not sure how to best handle this yet

## DIVISION ##
# not sure how to best handle this yet

#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::ExtendedStateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::ExtendedStateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(sys.D)
Base.size(sys::ExtendedStateSpace, d) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:ExtendedStateSpace} = S

function Base.getindex(sys::ExtendedStateSpace, inds...)
    if size(inds, 1) != 2
        error("Must specify 2 indices to index statespace model")
    end
    rows, cols = index2range(inds...) # FIXME: ControlSystems.index2range(inds...)
    return ExtendedStateSpace(
        copy(sys.A),
        sys.B[:, cols],
        sys.C[rows, :],
        sys.D[rows, cols],
        sys.Ts,
    )
end

#####################################################################
##                        Display Functions                        ##
#####################################################################

Base.print(io::IO, sys::ExtendedStateSpace) = show(io, sys)

function Base.show(io::IO, sys::ExtendedStateSpace)
    # Compose the name vectors
    println(io, typeof(sys))
    if nstates(sys) > 0
        println(io, "A = \n", _string_mat_with_headers(sys.A))
        println(io, "B1 = \n", _string_mat_with_headers(sys.B1))
        println(io, "B2 = \n", _string_mat_with_headers(sys.B2))
        println(io, "C1 = \n", _string_mat_with_headers(sys.C1))
        println(io, "C2 = \n", _string_mat_with_headers(sys.C2))
        println(io, "D11 = \n", _string_mat_with_headers(sys.D11))
        println(io, "D12 = \n", _string_mat_with_headers(sys.D12))
        println(io, "D21 = \n", _string_mat_with_headers(sys.D21))
        println(io, "D22 = \n", _string_mat_with_headers(sys.D22))
    else
        println(io, "The extended statespece model has no states..!")
    end

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
