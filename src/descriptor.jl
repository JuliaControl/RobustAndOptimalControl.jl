
"""
    DescriptorSystems.dss(sys::AbstractStateSpace)

Convert `sys` to a descriptor statespace system from [DescriptorSystems.jl](https://andreasvarga.github.io/DescriptorSystems.jl/dev/index.html)
"""
function DescriptorSystems.dss(sys::AbstractStateSpace)
    A,B,C,D = ssdata(sys)
    DescriptorSystems.dss(A, B, C, D; Ts = isdiscrete(sys) ? sys.Ts : 0)
end

function ControlSystemsBase.ss(G::DescriptorSystems.DescriptorStateSpace)
    try
        G,r = DescriptorSystems.gss2ss(G)
        if r < size(G.A, 1)
            G = DescriptorSystems.dss2ss(G, fast=false)
        end
    catch
        G,_ = DescriptorSystems.dss2ss(G, simple_infeigs = false, fast=false)
    end
    G.Ts > 0 ? ss(G.A, G.B, G.C, G.D, G.Ts) : ss(G.A, G.B, G.C, G.D)
end

"""
    n, ω = hinfnorm2(sys::LTISystem; kwargs...)

A numerically robust version of `hinfnorm` using DescriptorSystems.jl

For keyword arguments, see the docstring of `DescriptorSystems.ghinfnorm`, reproduced below
$(@doc(DescriptorSystems.ghinfnorm))
"""
function hinfnorm2(sys::LTISystem; kwargs...)
    sys, _ = balance_statespace(sys)
    DescriptorSystems.ghinfnorm(dss(ss(sys)); kwargs...)
end

function linfnorm2(sys::LTISystem; kwargs...)
    sys, _ = balance_statespace(sys)
    DescriptorSystems.glinfnorm(dss(ss(sys)); kwargs...)
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
    sysr, hs = baltrunc2(sys::LTISystem; residual=false, n=missing, kwargs...)

Compute the a balanced truncation of order `n` and the hankel singular values

For keyword arguments, see the docstring of `DescriptorSystems.gbalmr`, reproduced below
$(@doc(DescriptorSystems.gbalmr))
"""
function baltrunc2(sys::LTISystem; residual=false, n=missing, kwargs...)
    sysr, hs = DescriptorSystems.gbalmr(dss(sys); matchdc=residual, ord=n, kwargs...)
    ss(sysr), hs
end

"""
    sysr, hs, info = baltrunc_coprime(sys; residual = false, n = missing, factorization::F = DescriptorSystems.gnlcf, kwargs...)

Compute a balanced truncation of the left coprime factorization of `sys`.
See [`baltrunc2`](@ref) for additional keyword-argument help.

Coprime-factor reduction performs a coprime factorization of the model into \$P(s) = M(s)^{-1}N(s)\$ where \$M\$ and \$N\$ are stable factors even if \$P\$ contains unstable modes. After this, the system \$NM = \\begin{bmatrix}N & M \\end{bmatrix}\$ is reduced using balanced truncation and the final reduced-order model is formed as \$P_r(s) = M_r(s)^{-1}N_r(s)\$. For this method, the Hankel signular values of \$NM\$ are reported and the reported errors are \$||NM - N_rM_r||_\\infty\$. This method is of particular interest in closed-loop situations, where a model-reduction error \$||NM - N_rM_r||_\\infty\$ no greater than the normalized-coprime margin of the plant and the controller, guaratees that the closed loop remains stable when either \$P\$ or \$K\$ are reduced. The normalized-coprime margin can be computed with `ncfmargin(P, K)` ([`ncfmargin`](@ref)).

# Arguments:
- `factorization`: The function to perform the coprime factorization. A non-normalized factorization may be used by passing `RobustAndOptimalControl.DescriptorSystems.glcf`.
- `kwargs`: Are passed to `DescriptorSystems.gbalmr`, the docstring of which is reproduced below:
$(@doc(DescriptorSystems.gbalmr))
"""
function baltrunc_coprime(sys, info=nothing; residual=false, n=missing, factorization::F = DescriptorSystems.gnlcf, kwargs...) where F
    if info !== nothing && hasproperty(info, :NM)
        @unpack N, M, NM = info
    else
        N,M = factorization(dss(sys))
        A,E,B,C,D = DescriptorSystems.dssdata(N)
        NM = DescriptorSystems.dss(A,E,[B M.B],C,[D M.D])
    end
    NMr, hs = DescriptorSystems.gbalmr(NM; matchdc=residual, ord=n, kwargs...)
    
    A,E,B,C,D = DescriptorSystems.dssdata(DescriptorSystems.dss2ss(NMr)[1])
    
    nu = size(N.B, 2)
    BN = B[:, 1:nu]
    DN = D[:, 1:nu]
    BM = B[:, nu+1:end]
    DMi = pinv(D[:, nu+1:end])
    
    Ar = A - BM * (DMi * C)
    Cr = (DMi * C)
    Br = BN  - BM * (DMi * DN)
    Dr = (DMi * DN)

    ss(Ar,Br,Cr,Dr,sys.timeevol), hs, (; NM, N, M, NMr)
end


"""
    baltrunc_unstab(sys::LTISystem; residual = false, n = missing, kwargs...)

Balanced truncation for unstable models. An additive decomposition of sys into `sys = sys_stable + sys_unstable` is performed after which `sys_stable` is reduced. The order `n` must not be less than the number of unstable poles.

See `baltrunc2` for other keyword arguments.
"""
function baltrunc_unstab(sys::LTISystem, info=nothing; residual=false, n=missing, kwargs...)
    if info !== nothing && hasproperty(info, :stab)
        @unpack stab, unstab = info
    else
        stab, unstab = DescriptorSystems.gsdec(dss(sys); job="stable", kwargs...)
    end
    nx_unstab = size(unstab.A, 1)
    if n isa Integer && n < nx_unstab
        error("The model contains $(nx_unstab) poles outside the stability region, the reduced-order model must be of at least this order.")
    end
    sysr, hs = DescriptorSystems.gbalmr(stab; matchdc=residual, ord=n-nx_unstab, kwargs...)
    ss(sysr + unstab), hs, (; stab, unstab)
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