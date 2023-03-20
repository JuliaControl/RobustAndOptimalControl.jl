# When are two systems similar?
What does it mean for two systems to be *similar* to each other? Well, it depends on how we measure similarity. How we chose to measure similarity between two systems, or two models of the same system, may depend on what we are going to use the models for. In this little example, borrowed from the nice book
> Feedback Systems: An Introduction for Scientists and Engineers, by Åström and Murray

we will consider three different system models:
```math
P_1 = \dfrac{100}{s + 1}, \quad P_2 = \dfrac{100}{(s+1)(0.0025s+1)^2}, \quad P_3 = \dfrac{100}{s-1}
```

```@example SIMILARITY
using ControlSystemsBase, RobustAndOptimalControl, Plots
s = tf("s")
P1 = 100 / (s + 1)
P2 = 100 / ((s+1)*(0.025s+1)^2)
P3 = 100 / (s - 1)
P1, P2, P3
```

We start by having a look at the step responses of models $P_1$ and $P_2$:
```@example SIMILARITY
plot(step.([P1, P2], 6), lab=["P1" "P2"])
```
They sure look very *similar*, don't they? If we observed some noisy data from an experiment that look something like this
```@example SIMILARITY
res = step(P1, 0:0.1:6)
scatter!(res.t, res.y' .+ 2 .* randn.(), lab="Noisy observed data")
```
it would be very difficult to tell which of the two models $P_1$ and $P_2$ that best fit the data. However, if we now close the loop around these two "similar" models, we get the following:
```@example SIMILARITY
plot(step.([feedback(P1), feedback(P2)], 0.3))
```
Wow, that's a pretty big difference! Closing the loop around ``P_1`` resulted in a stable system, while closing the loop around ``P_2`` resulted in an unstable system! Surely, these two models cannot be considered similar to each other? 

Let's move on and compare the step responses of models $P_1$ and $P_3$:
```@example SIMILARITY
plot(step.([P1, P3], 2), lab=["P1" "P2"])
```
These two are obviously very *dissimilar*, one being stable and the other one not, yet, when we close the loop around these two models, we get the following:
```@example SIMILARITY
plot(step.([feedback(P1,1), feedback(P3,1)], 0.3))
```
we get very similar step responses! So, what does it *really mean* for two systems to be similar to each other? Simply looking at a simulation of the system might not always be sufficient. Let's have a look at the classical Bode and Nyquist curves:
```@example SIMILARITY
plot(
    bodeplot([P1, P2, P3], lab=["P1" "P2" "P3"]),
    nyquistplot([P1, P2, P3], lab=["P1" "P2" "P3"], xlims=(-3,1), ylims=(-3,1)),
)
```
Interestingly, all three models have the same gain for low frequencies, but the phase curves, and thus also the Nyquist curves, differ a lot. The Nyquist curve gives us an intuitive indication of how the system will perform when we close the loop. Here, the two models that are similar to each other are $P_1$ and $P_3$ (at least with these axis limits), while the model $P_2$ clearly encircle the critical point -1 in an unfortunate way.

What measure of similarity could we then use that takes into account how the system will perform when we close the loop? If we have a look at a standard similarity measure such as the ``H_\infty `` norm, we get that the models ``P_1`` and ``P_2`` are somewhat similar to each other, while ``P_1`` and ``P_3`` are not:
```@example SIMILARITY
hinfnorm(P1 - P2)[1], hinfnorm(P1 - P3)[1]
```
this does not align at all with how the systems behaved under feedback. Another metric, suitable for measuring similarity between systems when they are used in feedback, is the ``\nu``-gap metric:

```@example SIMILARITY
nugap(P1, P2)[1], nugap(P1, P3)[1]
```
this metric is always between 0 and 1, and a small "gap" indicates that the two models compared are similar. This metric aligns much better with how the systems behaved under feedback, indicating that ``P_1`` and ``P_3`` are similar to each other, while ``P_1`` and ``P_2`` are not.

This metric has an interesting relation to the normalized-coprime factor margin, [`ncfmargin`](@ref):
If controller ``K`` stabilizes ``P`` with an `ncfmargin(P, K)` ``= m``, then ``K`` will also stabilize all systems ``P'`` with a ``\nu``-gap from ``P`` of at most ``m``. This means that if our model error is small in the sense that ``\nu``-gap is small, and we design a controller with a large NCF-margin, then we can be confident that the controller will still stabilize the system even if the model is not perfect!

This property of the ``\nu``-gap metric and the NCF-margin is very useful in practice, and can be used for model-order reduction with guaranteed preservation of stability etc. see [`baltrunc_coprime`](@ref) and [Model-order reduction for uncertain models](@ref) for more info on this topic.

## Summary
This example has demonstrated that what it means for two models to be similar to each other might not always be a straightforward question to answer. Models that have very similar step responses, and simulation characteristics in general, might behave dramatically different when placed in a feedback loop. 
