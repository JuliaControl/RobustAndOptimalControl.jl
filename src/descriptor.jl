
"""
    DescriptorSystems.dss(sys::AbstractStateSpace)

Convert `sys` to a descriptor statespace system from [DescriptorSystems.jl](https://andreasvarga.github.io/DescriptorSystems.jl/dev/index.html)
"""
function DescriptorSystems.dss(sys::AbstractStateSpace)
    A,B,C,D = ssdata(sys)
    DescriptorSystems.dss(A, B, C, D; Ts = isdiscrete(sys) ? sys.Ts : 0)
end

function ControlSystems.ss(G::DescriptorSystems.DescriptorStateSpace)
    try
        G,r = DescriptorSystems.gss2ss(G)
        if r < size(G.A, 1)
            G = DescriptorSystems.dss2ss(G, fast=false)
        end
    catch
        G,_ = DescriptorSystems.dss2ss(G, simple_infeigs = false, fast=false)
    end
    G.Ts >= 0 ? ss(G.A, G.B, G.C, G.D) : ss(G.A, G.B, G.C, G.D, G.Ts)
end

"""
    n, ω = hinfnorm2(sys::LTISystem; kwargs...)

A numerically robust version of `hinfnorm` using DescriptorSystems.jl

For keyword arguments, see the docstring of `DescriptorSystems.ghinfnorm`, reproduced below
$(@doc(DescriptorSystems.ghinfnorm))
"""
function hinfnorm2(sys::LTISystem; kwargs...)
    DescriptorSystems.ghinfnorm(dss(ss(sys)); kwargs...)
end

"""
    n = h2norm(sys::LTISystem; kwargs...)

A numerically robust version of `norm` using DescriptorSystems.jl

For keyword arguments, see the docstring of `DescriptorSystems.gh2norm`, reproduced below
$(@doc(DescriptorSystems.gh2norm))
"""
function h2norm(sys::LTISystem; kwargs...)
    DescriptorSystems.gh2norm(dss(ss(sys)); kwargs...)
end

"""
    n, hsv = hankelnorm(sys::LTISystem; kwargs...)

Compute the hankelnorm and the hankel singular values

For keyword arguments, see the docstring of `DescriptorSystems.ghanorm`, reproduced below
$(@doc(DescriptorSystems.ghanorm))
"""
function hankelnorm(sys::LTISystem; kwargs...)
    DescriptorSystems.ghanorm(dss(ss(sys)); kwargs...)
end

"""
    nugap(sys0::LTISystem, sys1::LTISystem; kwargs...)

Compute the ν-gap metric between two systems. See also [`ncfmargin`](@ref).

For keyword arguments, see the docstring of `DescriptorSystems.gnugap`, reproduced below
$(@doc(DescriptorSystems.gnugap))
"""
function nugap(sys0::LTISystem, sys1::LTISystem; kwargs...)
    DescriptorSystems.gnugap(dss(ss(sys0)), dss(ss(sys1)); kwargs...)
end

const νgap = nugap



"""
    baltrunc2(sys::LTISystem; residual=false, n=missing, kwargs...)

Compute the a balanced truncation of order `n` and the hankel singular values

For keyword arguments, see the docstring of `DescriptorSystems.gbalmr`, reproduced below
$(@doc(DescriptorSystems.gbalmr))
"""
function baltrunc2(sys::LTISystem; residual=false, n=missing, kwargs...)
    sysr, hs = DescriptorSystems.gbalmr(dss(sys); matchdc=residual, ord=n, kwargs...)
    ss(sysr), hs
end

# function coprime_baltrunc(sys; residual=false, n=missing, kwargs...)
#     N,M = DescriptorSystems.glcf(dss(sys), mindeg=true, mininf=true, fast=false)
#     nu = size(N.B, 2)
#     # A,E,B,C,D = DescriptorSystems.dssdata(N)
#     # NM = DescriptorSystems.dss(A,E,[B M.B],C,[D M.D])
#     # Nr1, hs = DescriptorSystems.gbalmr(N; matchdc=residual, ord=n-size(M.A, 1), kwargs...)
#     NM = [N M]
#     sysr, hs = DescriptorSystems.gbalmr(NM; matchdc=residual, ord=n, kwargs...)

#     A,E,B,C,D = DescriptorSystems.dssdata(sysr)

#     Nr = DescriptorSystems.dss(A,E,B[:, 1:nu],C,D[:, 1:nu])
#     Mr = DescriptorSystems.dss(A,E,B[:, nu+1:end],C,D[:, nu+1:end])
#     sysr = Mr \ Nr
#     ss(sysr), hs
#     # Nr, Mr
# end


"""
    baltrunc_unstab(sys::LTISystem; residual = false, n = missing, kwargs...)

Balanced truncation for unstable models. An additive decomposition of sys into `sys = sys_stable + sys_unstable` is performed after which `sys_stable` is reduced. The order `n` must not be less than the number of unstable poles.

See `baltrunc2` for other keyword arguments.
"""
function baltrunc_unstab(sys::LTISystem; residual=false, n=missing, kwargs...)
    stab, unstab = DescriptorSystems.gsdec(dss(sys); job="stable", kwargs...)
    nx_unstab = size(unstab.A, 1)
    if n isa Integer && n < nx_unstab
        error("The model contains $(nx_unstab) poles outside the stability region, the reduced-order model must be of at least this order.")
    end
    sysr, hs = DescriptorSystems.gbalmr(stab; matchdc=residual, ord=n-nx_unstab, kwargs...)
    ss(sysr + unstab), hs
end

"""
    stab, unstab = stab_unstab(sys; kwargs...)

Decompose `sys` into `sys = stab + unstab` where `stab` contains all stable poles and `unstab` contains unstable poles. See $(@doc(DescriptorSystems.gsdec)) for keyword arguments (argument `job` is set to `"stable"` in this function).
"""
function stab_unstab(sys; kwargs...)
    stab, unstab = DescriptorSystems.gsdec(dss(sys); job="stable", kwargs...)
    ss(stab), ss(unstab)
end

##

function Base.:\(G1::AbstractStateSpace, G2::AbstractStateSpace)
    g1,g2 = dss(G1), dss(G2)
    try
        ss(g1\g2)
    catch
        error("Failed to solve G1\\G2, the product might not be proper. To keep the product on descriptor form, use dss(G1)\\dss(G2)")
    end
end