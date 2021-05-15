import ControlSystems as CS
import ControlSystems: nstates, blockdiag
function named_ss(sys;
    state_names = [Symbol("x$i") for i in 1:sys.nx],
    input_names = [Symbol("u$i") for i in 1:sys.nu],
    output_names = [Symbol("y$i") for i in 1:sys.ny],
    )
    length(state_names) == sys.nx  || throw(ArgumentError("Length of state names must match sys.nx"))
    length(input_names) == sys.nx  || throw(ArgumentError("Length of input names must match sys.nx"))
    length(output_names) == sys.ny || throw(ArgumentError("Length of state names must match sys.ny"))
    A,B,C,D = ssdata(sys)

    # ax_x = ViewAxis(1:sys.nx, (; (state_names .=> (1:sys.nx))...))
    ax_x = Axis(state_names...)
    ax_u = Axis(input_names...)
    ax_y = Axis(output_names...)

    An = ComponentMatrix(A, ax_x, ax_x)
    Bn = ComponentMatrix(B, ax_x, ax_u)
    Cn = ComponentMatrix(C, ax_y, ax_x)
    Dn = ComponentMatrix(D, ax_y, ax_u)

    HeteroStateSpace(An, Bn, Cn, Dn, sys.timeevol)
end


# function ComponentArrays.Axis(names::Symbol...)
#     Axis((; (names .=> (1:length(names)))...))
# end

function ComponentArrays.Axis(axs::Axis...)
    names_ = [[keys(ComponentArrays.indexmap(ax))...] for ax in axs]
    names = reduce(vcat, names_)
    length(unique(names)) == length(names) || throw(ArgumentError("Axis names repeated: $names"))
    Axis(names...)
end

function Base.hcat(x1::ComponentMatrix,x2::ComponentMatrix,x3::ComponentMatrix,x::ComponentMatrix...) 
    x = [x1,x2,x3,x...]
    X = getdata.(x)
    axs = getaxes.(x)
    ax1s = first.(axs)
    ax2s = last.(axs)
    # all(==(ax1s[1]), ax1s) || throw(ArgumentError("Cannot hcat of not all row-axes are the same. $ax1s"))
    ax_1 = getaxes(x[1])[1]
    ax_2 = Axis(ax2s...)
    ComponentMatrix(reduce(hcat, X), ax_1, ax_2)
end

function Base.vcat(x1::ComponentMatrix,x2::ComponentMatrix,x3::ComponentMatrix,x::ComponentMatrix...)
    x = [x1,x2,x3,x...]
    X = getdata.(x)
    axs = getaxes.(x)
    ax1s = first.(axs)
    ax2s = last.(axs)
    # all(==(ax2s[1]), ax2s) || throw(ArgumentError("Cannot vcat of not all column-axes are the same. $ax2s"))
    ax_2 = ax2s[1]
    ax_1 = Axis(ax1s...)
    ComponentMatrix(reduce(vcat, X), ax_1, ax_2)
end


function Base.:+(s1::HeteroStateSpace{TE, MT1}, s2::HeteroStateSpace{TE, MT2}) where {TE <: CS.TimeEvolution, MT1 <: ComponentMatrix, MT2 <: ComponentMatrix}
    #Ensure systems have same dimensions
    if size(s1) != size(s2)
        error("Systems have different shapes.")
    end
    timeevol = CS.common_timeevol(s1,s2)
    ax1 = getaxes(s1.A)[1]
    ax2 = getaxes(s2.A)[1]
    T = promote_type(eltype(s1.A),eltype(s2.A))
    A = [[s1.A                   ComponentMatrix(zeros(T, nstates(s1), nstates(s2)), ax1, ax2)];
    [ComponentMatrix(zeros(T, nstates(s2), nstates(s1)), ax2, ax1)    s2.A]]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C]
    D = s1.D + s2.D

    return HeteroStateSpace(A, B, C, D, timeevol)
end

function Base.:*(s1::HeteroStateSpace{TE, MT1}, s2::HeteroStateSpace{TE, MT2}) where {TE <: CS.TimeEvolution, MT1 <: ComponentMatrix, MT2 <: ComponentMatrix}
    #Check dimension alignment
    #Note: s1*s2 = y <- s1 <- s2 <- u
    if s1.nu != s2.ny
        error("s1*s2: s1 must have same number of inputs as s2 has outputs")
    end
    timeevol = CS.common_timeevol(s1,s2)
    ax1 = getaxes(s1.A)[1]
    ax2 = getaxes(s2.A)[1]
    axu1, axu2 = getaxes(s1.B)[2], getaxes(s2.B)[2]
    T = promote_type(eltype(s1.A),eltype(s2.A))
    A12 = s1.B*s2.C
    B1 = s1.B*s2.D
    axa12 = getaxes(B1)[1]
    A = [[s1.A    A12];
         [ComponentMatrix(zeros(T, s2.nx, s1.nx), ax2, ax1)  s2.A]]
    B = ComponentMatrix(getdata([B1; s2.B]), Axis(ax1, ax2), axu2)
    C = [s1.C   s1.D*s2.C]
    # D = ComponentMatrix(getdata(s1.D*s2.D), getaxes(s1.D)[1], getaxes(s2.D)[2])
    D = s1.D*s2.D

    return HeteroStateSpace(A, B, C, D, timeevol)
end

function Base.hcat(systems::HeteroStateSpace{<:TE, <:ComponentMatrix}...) where TE
    ny = systems[1].ny
    if !all(s.ny == ny for s in systems)
        error("All systems must have same output dimension")
    end
    A = blockdiag([s.A for s in systems])
    B = blockdiag([s.B for s in systems])
    C = hcat([s.C for s in systems]...)
    D = hcat([s.D for s in systems]...)
    timeevol = common_timeevol(systems...)
    return HeteroStateSpace(A, B, C, D, timeevol)
end

function Base.vcat(systems::HeteroStateSpace{<:TE, <:ComponentMatrix}...) where TE
    # Perform checks
    nu = systems[1].nu
    if !all(s.nu == nu for s in systems)
        error("All systems must have same input dimension")
    end
    A = blockdiag([s.A for s in systems])
    B = vcat([s.B for s in systems]...)
    C = blockdiag([s.C for s in systems])
    D = vcat([s.D for s in systems]...)
    timeevol = common_timeevol(systems...)
    return HeteroStateSpace(A, B, C, D, timeevol)
end

function ControlSystems.blockdiag(mats::Vector{<:ComponentMatrix{T}}) where T
    rows = Int[size(m, 1) for m in mats]
    cols = Int[size(m, 2) for m in mats]

    axs = getaxes.(mats)
    ax1 = Axis(first.(axs)...)
    ax2 = Axis(last.(axs)...)
    res = ComponentMatrix(zeros(T, sum(rows), sum(cols)), ax1, ax2)
    m = 1
    n = 1
    for ind=1:length(mats)
        mat = mats[ind]
        i = rows[ind]
        j = cols[ind]
        res[m:m + i - 1, n:n + j - 1] = mat
        m += i
        n += j
    end
    return res
end