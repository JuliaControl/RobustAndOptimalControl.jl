"""
    nugap(G; map = map)

Compute the νgap between the nominal system `Gₙ` represented by the first particle index in `G`, and all other systems in `G`.
Returns a `Particles` object with the νgap for each system in `G`.

See `with_nominal` to endow uncertain values with a nominal value, and `nominal` to extract the nominal value.

The value returned by this function, `νᵧ` is useful for robust synthesis, by designing a controller for the nominal system `Gₙ`, that achieves an [`ncfmargin`](@ref) of at least `νᵧ` is guaranteed to stabilize all realizations within `G`. 

To speed up computation for large systems, a threaded or distributed `map` function can be supplied, e.g., `ThreadTools.tmap` or `Distributed.pmap`.
"""
function RobustAndOptimalControl.nugap(G, i = 1; map = map)
    ControlSystemsBase.numeric_type(G) <: MonteCarloMeasurements.AbstractParticles ||
        throw(ArgumentError("When calling nugap with a single argument, G must be an uncertain system with particles where the nominal value is the first particle."))
    Gn = RobustAndOptimalControl.sys_from_particles(G, i)
    N = nparticles(G.A)
    gaps = map(1:N) do i
        Gi = RobustAndOptimalControl.sys_from_particles(G, i)
        nugap(Gn, Gi)[1]
    end
    Particles(gaps)
end



##

"""
    nu_reduction(G, g=0.1; gap = nugap(G))

Reduce the number of particles in an uncertain system `G` by removing all particles that are within the νgap `g` of the nominal system `Gₙ`.

Note: If `G` has a stochastic interpretation, i.e., the coefficients come from some distribution, this interpretation will be lost after reduction, mean values and standard deviations will not be preserved. The reduced system should instead be interpreted as preserving worst-case uncertainty.

If the `gap = nugap(G)` has already been precomputed, it can be supplied as an argument to avoid potentially costly recomputaiton.
"""
function nu_reduction(G, g=0.1; gap = nugap(G))
    nind = argmin(gap.particles)
    keepinds = [nind; findall(gap.particles .> g)] # keep also nominal
    syss = [RobustAndOptimalControl.sys_from_particles(G, i) for i in keepinds]
    RobustAndOptimalControl.ss2particles(syss)
end


"""
    nu_reduction_recursive(G, g = 0.1; gap = nugap(G), keepinds = Set{Int}(1), verbose = false)

Find a νgap cover of balls of radius `g` (in the νgap metric) that contain all realizations in `G`. 

If the `gap = nugap(G)` has already been precomputed, it can be supplied as an argument to avoid potentially costly recomputaiton. If a manually computed `gap` is supplied, you must also supply `keepinds=Set{Int}(index)` where `index` is the index of the nominal system in `G` used to compute `gap`.

The returned cover `Gr` is of the same type as `G`, but with a smaller number of particles. A controller designed for `Gr` that achieves a [`ncfmargin`](@ref) of at least `g` for all realizations in `Gr` will stabilize all realizations in the original `G`. The extreme case cover where `Gr = Gnominal` is a single realization only can be computed by calling `g = nugap(G, i)` where `i` is the index of the nominal system in `G`.

# Arguments:
- `G`: An uncertain model in the form of a `StateSpace{TE, Particles}` (a multi-model).
- `g`: The radius of the balls in the νgap cover.
- `gap`: An optional precomputed gap
- `verbose`: Print progress
"""
function nu_reduction_recursive(G, g=0.1; gap = nugap(G), discardinds=Set{Int}(), keepinds=Set{Int}(1), verbose = false)
    nparticles(gap) >= 10 && g <= 0.01 && @warn "This algorithm has a complexity of O(N²) for small `g`" maxlog=1
    discardinds = union(discardinds, Set(findall(gap.particles .< g)))
    verbose && println("Covered $(length(discardinds)) particles")
    maxind, _ = argmax(enumerate(gap.particles)) do (i,g)
        if i ∈ discardinds
            zero(g)
        else
            g
        end
    end
    # @show maxind
    if maxind ∈ discardinds
        N = nparticles(G.A)
        keepers = sort(unique([collect(keepinds); [i for i in 1:N if i ∉ discardinds]]))
        syss = [RobustAndOptimalControl.sys_from_particles(G, i) for i in keepers]
        return RobustAndOptimalControl.ss2particles(syss)
    else
        push!(keepinds, maxind)
    end
    gap = nugap(G, maxind)
    nu_reduction_recursive(G, g; gap, discardinds, keepinds, verbose)
end

