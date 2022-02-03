# Extended statespace systems

An [`ExtendedStateSpace`](@ref) system represents a two-input, two-output system.

A [`ChannelStateSpace`](@ref) represents a statespace system with open and closed input-output channels, like the example shown in the block diagram below. There is one open channel, $w \rightarrow z$, and 4 closed channels. Each closed channel is closed through an [`AbstractBlock`](@ref). The blocks in the diagram are just examples. The `ChannelStateSpace` type represents the entire structure in the block diagram, including $P$ and all the blocks. 
```
z         ┌─────────┐          w
◄─────────┤         │◄────────── # Performance channel (open)
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
│ │ │ └─────►│ K ├───────┘ │ │ │ # Feedback channel
│ │ │        └───┘         │ │ │
│ │ │                      │ │ │
│ │ │        ┌───┐         │ │ │
│ │ └───────►│ Δ ├─────────┘ │ │ # Uncertainty channel
│ │          └───┘           │ │
│ │                          │ │
│ │          ┌───┐           │ │
│ └─────────►│ τ ├───────────┘ │ # Delay channel
│            └───┘             │
│                              │
│            ┌───┐             │
└───────────►│LPV├─────────────┘ # Linear parameter-varying channel
             └───┘
```
Internally, [`ChannelStateSpace`](@ref) maintains a vector of [`IOChannel`](@ref)s.

## Channels
An [`IOChannel`](@ref) represents a mapping from inputs to outputs through a dynamical system $P$. A channel has the following properties:
- A type parameter `T` in `IOChannel{T}` that specifies the type of the channel.
- A name (`Symbol`) for the entire group of inputs `u`, and one name for the group of outputs `y`.
- A vector of the names of all individual inputs `us` and one vector for the outputs `ys`.
- A vector of `blocks` through which the channel is closed. If the channel is open, this vector is empty.

Use [`isopen`](@ref), [`isclosed`](@ref) to figure out if a channel is open or closed. 

Channels are always sorted internally such that channels that are open are stored first, and closed channels last (technicalities applies, see below). In this way, the `ChannelStateSpace` type can be seen as a specialization of `ExtendedStateSpace` where there is an additional partitioning of the two input and output channels. 

The following channels are defined in this package
```@example 
subtypes(AbstractChannelType)
```

A channel *does not* include the connection and dynamics matrices, i.e., the tuple $A,B,C,D$, since there are cross terms between channels, like $D_{12}$ in the [`ExtendedStateSpace`](@ref) type, that could not be represented like this. Instead, the [`ChannelStateSpace`](@ref) maintains an inner `ExtendedStateSpace` object, and the channels are placed in the upper or lower channel of the `ExtendedStateSpace` based on their algebraic properties, more on this in [Algebraic behavior of channels](@ref)

## Blocks
A closed channel is closed through an `AbstractBlock`. Typical examples include
- Uncertainty channel, closed through an [`UncertainElement`](@ref).
- Delay channel, closed through pure delay elements (no type for delay element exists yet). 


## Algebraic behavior of channels
When two systems `s1,s2` are connected in series `s1*s2`, standard statespace systems connect the outputs of `s2` to the inputs of `s1`, and the resulting product system has the outputs of `s1` and the inputs of `s2`. In the `ChannelStateSpace` type, this type of behavior is represented by the `DefaultChannel`. This is an open channel that behaves just like an ordinary `StateSpace` system.

Closed-channels, on the other hand, do not behave like this, rather, if two systems with closed channels of the same type are connected in series, the closed channel grows in size to include the input-output mappings of *both* systems, and the corresponding channel of the product system contains the blocks from both systems. 

*Typically*, open channels behave like the standard `StateSpace` type, and closed channels "append" rather than series, but there are some exceptions:
- A [`ConstrainedChannel`](@ref) represents outputs that a constrained to lie within some constraint set, e.g., for MPC applications. A [`ConstrainedChannel`](@ref) is not closed around anything, and is thus to be considered open, but when two systems are connected, all constrained outputs should be outputs of the connected system. 
An example of when this situation comes up is given below.

The system $P$ contains measured outputs and constrained outputs, $y$ and $v$. Some of the inputs $u$ of $P$ are directly fed through to $v$ to indicate input constraints. Loop shaping on $P$ with pre and post-compensators $W_1, W_2$, forms the system
$$P_s = W_2 P W_1$$
which, if $P$ was a standard `StateSpace`, would have new inputs and outputs $u_s, y_s$. However, if we have specified constraints on the output $v$ of $P$, these will in general not hold for the outputs of the scaled plant $P_s$. Hence, $P_s$ must include the original constrained outputs $v$ among it's outputs.

Since the open/closedness of a channel doesn't always indicate the algebraic properties of a channel, we instead make use of a trait-based solution. Each channel type defines an implementation of the function [`algebraic_trait`](@ref) which returns either `SeriesTrait` or `AppendTrait`, indicating its behavior when systems are multiplied. The inner [`ExtendedStateSpace`](@ref) of a `ChannelStateSpace` is thus technically not partitioned based on open and closed channels, rather, it's partitioned based on the `algebraic_trait`.

### Promotion
Before two systems can interact through algebraic operations, they must be promoted to a common supertype. In this case, the supertype is the `ChannelStateSpace` that contains the union of all the channel types of both systems. The channels that were added in the promotion step are simply empty, i.e., contains no signals and no blocks.

`LTISystem` objects that are not instances of `ChannelStateSpace` are promoted to a `ChannelStateSpace` with the `DefaultChannel` type. 


### LFT on ChannelStateSpace
The literature commonly talks about "upper" and "lower" linear fractional transforms. The terminology upper/lower comes from the TITO system view (`ExtendedStateSpace`) where the lower LFT (`lft(P, K, :l)`) closes the lower loop around `K`. For `ChannelStateSpace` types, the terms upper and lower have no meaning since we do no longer have the TITO convention, instead, channels are typed and we refer to the type of the channel we want to close in the LFT, e.g.,
```julia
lft(P, K, DefaultChannel)
```
closes the default channel over $K$. In this scenario, the `DefaultChannel` is a bit special, since `K` is a regular `LTISystem`, the resulting "closed-loop" system no longer has a `DefaultChannel` (we never keep standard `LTISystem`s as blocks of a `ChannelStateSpace`). You can only call `lft` on open channels, since closed channels are by definition already a closed LFT, `lft(P, block, ch)`.


## Creating ChannelStateSpace
Any `LTISystem` can be converted to a `ChannelStateSpace` by calling the constructor shorthand `css(sys)`, this gives them a `DefaultChannel`. Special kinds of `ChannelStateSpace` systems have their own constructors, such as an uncertain statespace model, with the constructor `uss`. 

If you would like to manually create a `ChannelStateSpace` system with specified channels, call
```julia
css(sys, channels...)
```
where `channels::IOChannel...` contain indices into the system `sys`. Indices that are not present in any of the `channels` will be placed in the default channel.