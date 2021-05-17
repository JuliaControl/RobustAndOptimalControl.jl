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
            error("D21 must have the same column size as B1")
        elseif nu != size(D12, 2)
            error("D12 must have the same column size as B2")
        elseif nu != size(D22, 2)
            error("D22 must have the same column size as B2")
        elseif nz != size(D11, 1)
            error("D11 must have the same row size as C1")
        elseif nz != size(D12, 1)
            error("D12 must have the same row size as C1")
        elseif ny != size(D21, 1)
            error("D21 must have the same row size as C2")
        elseif ny != size(D22, 1)
            error("D22 must have the same row size as C2")
        end


        new{TE,T}(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
    end
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
    D11::AbstractArray = zeros(T, size(C1, 1), size(B1, 2)),
    D12::AbstractArray = zeros(T, size(C1, 1), size(B2, 2)),
    D21::AbstractArray = zeros(T, size(C2, 1), size(B1, 2)),
    D22::AbstractArray = zeros(T, size(C2, 1), size(B2, 2)),
    Ts = nothing,
) where T
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

ControlSystems.StateSpace(s::ExtendedStateSpace) = ss(ssdata(s)..., s.timeevol)

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

## NEGATION ##
function Base.:-(sys::ST) where ST <: ExtendedStateSpace
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(sys)
    ST(A, B1, B2, -C1, -C2, D11, D12, D21, D22, sys.timeevol)
end

#####################################################################
##                       Indexing Functions                        ##
#####################################################################
Base.ndims(::ExtendedStateSpace) = 2 # NOTE: Also for SISO systems?
Base.size(sys::ExtendedStateSpace) = (noutputs(sys), ninputs(sys)) # NOTE: or just size(sys.D)
Base.size(sys::ExtendedStateSpace, d) = d <= 2 ? size(sys)[d] : 1
Base.eltype(::Type{S}) where {S<:ExtendedStateSpace} = S
ControlSystems.numeric_type(sys::ExtendedStateSpace) = eltype(sys.A)

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

function ControlSystems.lft(G::ExtendedStateSpace, Δ, type=:l)

    # if !(G.nu > Δ.ny && G.ny > Δ.nu)
    #     error("Must have G.nu > Δ.ny and G.ny > Δ.nu for lower/upper lft")
    # end

    if type === :l
        feedback(G, Δ)
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

function ControlSystems.feedback(s1::ExtendedStateSpace, s2::ExtendedStateSpace;
    Wperm=:, Zperm=:)

    timeevol = common_timeevol(s1,s2)
    s1_B1 = s1.B1
    s1_B2 = s1.B2
    s1_C1 = s1.C1
    s1_C2 = s1.C2
    s1_D11 = s1.D11
    s1_D12 = s1.D12
    s1_D21 = s1.D21
    s1_D22 = s1.D22


    # s2_B1 = s2.B2 # These are reversed
    # s2_B2 = s2.B1
    # s2_C1 = s2.C2
    # s2_C2 = s2.C1
    # s2_D11 = s2.D22
    # s2_D12 = s2.D21
    # s2_D21 = s2.D12
    # s2_D22 = s2.D11

    s2_B1 = s2.B1 
    s2_B2 = s2.B2
    s2_C1 = s2.C1
    s2_C2 = s2.C2
    s2_D11 = s2.D11
    s2_D12 = s2.D12
    s2_D21 = s2.D21
    s2_D22 = s2.D22

    @sizecompat s1_B2 s2_D22
    @sizecompat s1_B2 s2_C2
    @sizecompat s2_D22 s1_D21
    @sizecompat s2_D22 s1_C2
    @sizecompat s1_D12 s2_C2
    @sizecompat s2_D12 s1_D22
    @sizecompat s1_D22 s2_C2

    if iszero(s1_D22) || iszero(s2_D22)
        A = [s1.A + s1_B2*s2_D22*s1_C2        s1_B2*s2_C2;
                s2_B2*s1_C2            s2.A + s2_B2*s1_D22*s2_C2]

        B1 = [
            s1_B1 + s1_B2*s2_D22*s1_D21
                    s2_B2*s1_D21            
        ]
        B2  = [
            s1_B2*s2_D21
            s2_B1 + s2_B2*s1_D22*s2_D21
        ]
        C1 = [s1_C1+s1_D12*s2_D22*s1_C2        s1_D12*s2_C2]
        C2 = [s2_D12*s1_C2           s2_C1+s2_D12*s1_D22*s2_C2]
        D11 = s1_D11 + s1_D12*s2_D22*s1_D21 
        D12 = s1_D12*s2_D21
        D21 = s2_D12*s1_D21           
        D22 = s2_D11 + s2_D12*s1_D22*s2_D21
    else
        R1 = try
            inv(I - s2_D22*s1_D22)
        catch
            error("Ill-posed feedback interconnection,  I - s2_D22*s1_D22 or I - s2_D22*s1_D22 not invertible")
        end

        R2 = try
            inv(I - s1_D22*s2_D22)
        catch
            error("Ill-posed feedback interconnection,  I - s2_D22*s1_D22 or I - s2_D22*s1_D22 not invertible")
        end

        A = [s1.A + s1_B2*R1*s2_D22*s1_C2        s1_B2*R1*s2_C2;
                s2_B2*R2*s1_C2            s2.A + s2_B2*R2*s1_D22*s2_C2]

        B1 = [s1_B1+s1_B2*R1*s2_D22*s1_D21;
                    s2_B2*R2*s1_D21]
        B2 = [s1_B2*R1*s2_D21; s2_B1 + s2_B2*R2*s1_D22*s2_D21]
        C1 = [s1_C1 + s1_D12*R1*s2_D22*s1_C2        s1_D12*R1*s2_C2]
        C2 = [s2_D12*R2*s1_C2           s2_C1+s2_D12*R2*s1_D22*s2_C2]
        D11 = [s1_D11 + s1_D12*R1*s2_D22*s1_D21]
        D12 = s1_D12*R1*s2_D21
        D21 = s2_D12*R2*s1_D21
        D22 = s2_D11 + s2_D12*R2*s1_D22*s2_D21
    end

    return ExtendedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22, timeevol)
end

ExtendedStateSpace(s::AbstractStateSpace) = ss(s.A, I(s.nx), s.B, I(s.nx), s.C; D22=s.D)

function system_mapping(P::ExtendedStateSpace)
    ss(P.A, P.B2, P.C2, P.D22)
end

function performance_mapping(P::ExtendedStateSpace)
    ss(P.A, P.B11, P.C11, P.D11)
end