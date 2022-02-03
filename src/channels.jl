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

const NamedIndex = Union{Symbol, Vector{Symbol}, Colon}

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

iterable(s::Symbol) = [s]
iterable(v) = v

names2indices(::Colon, allnames) = 1:length(allnames) 

function names2indices(names, allnames)
    inds = Union{Nothing, Int}[findfirst(==(n), allnames) for n in names]
    snames = string.(allnames)
    for i in eachindex(inds)
        if inds[i] === nothing
            # try finding symbols with given prefix
            newi = findall(startswith(string(names[i])), snames)
            newi === nothing || isempty(newi) && error("The indexed NamedSystem has no signal named $(names[i]), available names are $(allnames)")
            deleteat!(inds, i)
            foreach(j->insert!(inds, j+i-1, newi[j]), 1:length(newi))
        end
    end
    inds
end
function names2indices(name::Symbol, allnames)
    i = findfirst(==(name), allnames)
    if i === nothing
        # try finding symbols with given prefix
        i = findall(startswith(name), allnames)
        error("The indexed NamedSystem has no signal named $name, available names are $(allnames)")
    else
        i:i # return a vector rather than scalar for slices of matrices to not drop dim
    end
end

"""
A system P with multiple IOChannels closed around AbstractBlocks.
"""
channel_block2 = """
z         ┌─────────┐          w
◄─────────┤         │◄──────────
          │         │
┌─────────┤         │◄─────────┐
│         │         │          │
│ ┌───────┤    P    │◄───────┐ │
│ │       │         │        │ │
│ │ ┌─────┤         │◄─────┐ │ │
│ │ │     │         │      │ │ │
│ │ │y┌───┤         │◄───┐u│ │ │
│ │ │ │   └─────────┘    │ │ │ │
│ │ │ │                  │ │ │ │
│ │ │ │      ┌───┐       │ │ │ │
│ │ │ └─────►│ K ├───────┘ │ │ │
│ │ │        └───┘         │ │ │
│ │ │                      │ │ │
│ │ │        ┌───┐         │ │ │
│ │ └───────►│ Δ ├─────────┘ │ │
│ │          └───┘           │ │
│ │                          │ │
│ │          ┌───┐           │ │
│ └─────────►│ τ ├───────────┘ │
│            └───┘             │
│                              │
│            ┌───┐             │
└───────────►│LPV├─────────────┘
             └───┘
"""


abstract type AbstractBlock end

# QUESTION: is a block a dynamical system? I think it is since it's included in the feedback loop, but do we want the behavior of LTISystem?

# ==============================================================================
## Channel types ===============================================================
# ==============================================================================
abstract type AbstractChannelType end

struct DefaultChannel <: AbstractChannelType end
Base.@kwdef struct PerformanceChannel <: AbstractChannelType
    closed::Bool = false # This channel can be closed around a virtual performance block
end
struct LPVChannel <: AbstractChannelType end
struct UncertaintyChannel <: AbstractChannelType end
struct ConstrainedChannel <: AbstractChannelType end
struct DelayChannel <: AbstractChannelType end
struct AnalysisChannel <: AbstractChannelType end

# ==============================================================================
## Channel implementations =====================================================
# ==============================================================================
abstract type AbstractIOChannel end

const OPEN_TYPES = Union{DefaultChannel, ConstrainedChannel}
const CLOSED_TYPES = Union{LPVChannel, UncertaintyChannel, DelayChannel}
isopen(ch::OPEN_TYPES) = true
isclosed(ch::OPEN_TYPES) = false
isopen(ch::CLOSED_TYPES) = false
isclosed(ch::CLOSED_TYPES) = true


struct IOChannel{T} <: AbstractIOChannel
    u::Symbol
    y::Symbol
    us::Vector{Symbol}
    ys::Vector{Symbol}
    blocks::Vector{AbstractBlock}
end

# QUESTION: a block may have a size, so that one can not naively compare the lengths of u/y and the length of the blocks vector. Do we need an abstract block type?
# QUESTION: Should blocks be allowed to be complete ChannelStateSpace systems? then we could represent stuff like feedback lft(P, K) where K = lft(K0, lpv) etc., or should this kind of system be normalized into single LFT form immediately? Probably normalized since lft(A, lft(B, C)) is an lft(D, C), see skogestad appendix A.7.1
# QUESTION: should u/y be allowed to be anything or just vectors of symbolics? Could they be SignalGroups, or that is a separate thing? Probably the latter, keep types simple.
# TODO: Implements getproperty on IOCHannel so that we can be nu, ny, etc.


# ==============================================================================
## Channel traits ==============================================================
# ==============================================================================

## Multiplication trait ========================================================
abstract type AlgebraicTrait end
struct AppendTrait <: AlgebraicTrait end # All channels that are closed around blocks have the AppendTrait, like uncertainty channel
struct SeriesTrait <: AlgebraicTrait end # Open channels have the SeriesTrait, like Default channel
struct NoInteractionTrait <: AlgebraicTrait end # Channels that don't bother with each other, like Uncertainty and Constrained

algebraic_trait(c1::AbstractIOChannel{T1}, c2::AbstractIOChannel{T2}) where {T1 <: AbstractChannelType, T2 <: AbstractChannelType} = algebraic_trait(T1, T2)

algebraic_trait(::Type{T1}, ::Type{T2}) where {T1 <: AbstractChannelType, T2 <: AbstractChannelType} = NoInteractionTrait()

function Base.:*(s1::ChannelStateSpace, s2::ChannelStateSpace)
    for ch1 in s1.channels
        for ch2 in s2.channels
            mt = algebraic_trait(c1, ch2)
            channel_mul(mt, c1, c2, s1, s2)
        end
    end
end


#=
NOTE: For USS, the multiplication just looks like this
sys = invert_mappings(invert_mappings(s1)*invert_mappings(s2))
UncertainSS(sys, [s1.Δ; s2.Δ])
There are two things going on there, the invert_mapping puts the Δ block below,
the multiplication then series the top channel and appends the lower channel,
i.e., the Append action on the blocks is handled separately from the multiplication.
The two channels here are DefaultChannel and UncertaintyChannel. During the multiplication, standard series formulas are used for the DefaultChannel, i.e., the resulting A matrix is influenced by B1/C1 only. The resulting B matrix
will have only DefaultChannel influence in the B1 part, and the B2 part will become larger due to the append (B2 is not only append)
The same for C, the new C1 contains only C1/D11 while C2 is influenced by append

When multiplying many channels, one could thus sort/group the channels to create this exact TITO representation where one channel is series and one append, and then repartition the resulting system!
A good idea is probably to always keep the system in sorted group order based on the Channel multiplication trait such that the sorting is avoided. This would in a way maintain an internal TITO system, but where each input output channel in the TITO view is subdivided further

Keep all seres channels on top, and append/closed channels in the bottom. This aligns with the diagram and with * for ExtendedStateSpace
Maybe we can just keep an extended statespace in the ChannelStateSpace?

The ConstrainedChannel is technically open, but should behave like a closed channel in the sense that post multiplication W2*P should keep constraiend outputs of P as outputs of W2*P, hence, it's not obvious how to connect channels based on the open/closed properties alone. Perhaps this is just a special case of multiplying two systems where not all channels are present in the interface between them? What happens to the channels that are not connected? Are they discarded or appended to the inputs/outputs of the product? Probably the latter. Will this ever happen in practice? The type example is W2*P*W1 where P has contrained channel but W1, W2 only default channel
=#


## Addition trait ==============================================================


## LFT =========================================================================

function ControlSystems.lft(s1::ChannelStateSpace, block, channel)
    # We no longer have upper/lower lft, we instead specify the channel to close through block
end