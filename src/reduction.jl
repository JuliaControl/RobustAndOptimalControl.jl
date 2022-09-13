using ControlSystemsBase: ssdata

"""
    sysr, hs = frequency_weighted_reduction(G, Wo, Wi; residual=true)

Find Gr such that ||Wₒ(G-Gr)Wᵢ||∞ is minimized.
For a realtive reduction, set Wo = inv(G) and Wi = I.

If `residual = true`, matched static gain is achieved through "residualization", i.e., setting
```math
0 = A_{21}x_{1} + A_{22}x_{2} + B_{2}u
```
where indices 1/2 correspond to the remaining/truncated states respectively. This choice typically results in a better match in the low-frequency region and a smaller overall error.

Ref: Andras Varga and Brian D.O. Anderson, "Accuracy enhancing methods for the frequency-weighted balancing related model reduction"
https://elib.dlr.de/11746/1/varga_cdc01p2.pdf

Note: This function only handles exponentially stable models. To reduce unstable  and marginally stable models, decompose the system into stable and unstable parts using [`stab_unstab`](@ref), reduce the stable part and then add the unstable part back.
"""
function frequency_weighted_reduction(G, Wo, Wi, r=nothing; residual=true, atol=sqrt(eps()), rtol=1e-3)
    iscontinuous(G) || error("Discrete systems not supported yet.")
    A,B,C,D = ssdata(G)
    n = G.nx
    if Wo == 1 || Wo == I
        R = grampd(G, :o)
    else
        Wo = ss(Wo)
        if issiso(Wo) && size(G, 1) > 1
            Wo = Wo*I(size(G, 1))
        end
        Ao,Bo,Co,Do = ssdata(Wo)
        As = [
            A                                  zeros(size(A,1), size(Ao,2)) 
            Bo*C                               Ao                             
        ]
        Bs = [
            B
            Bo*D
        ]
        Cs = [
            Do*C   Co 
        ]
        Ds = Do*D
        WoG = ss(As, Bs, Cs, Ds)
        R = grampd(WoG, :o)[1:n, 1:n]
    end

    if Wi == 1 || Wi == I
        S = grampd(G, :c)
    else
        Wi = ss(Wi)
        if issiso(Wi) && size(G, 2) > 1
            Wi = Wi*I(size(G, 1))
        end
        Ai,Bi,Ci,Di = ssdata(ss(Wi))
        As = [
            A                                B*Ci
            zeros(size(Ai,1), size(A,2))     Ai
        ]
        Bs = [
            B*Di
            Bi
        ]
        Cs = [
            C   D*Ci
        ]
        Ds = D*Di
        GWi = ss(As, Bs, Cs, Ds)
        S = grampd(GWi, :c)[1:n, 1:n]
    end

    U,Σ,V = svd!(R*S)
    rmin = count(Σ .> sqrt(eps())*Σ[1])
    if r === nothing
        r = count(Σ .> max(atol,rtol*Σ[1]))
    end
    r = min(r, rmin)
    i1 = 1:r
    i2 = r+1:rmin
    U1 = U[:,i1]
    V1 = V[:,i1]
    U2 = U[:,i2]
    V2 = V[:,i2]
    Σ1 = Σ[i1]
    Y = Matrix(qr!(R'U1).Q)
    X = Matrix(qr!(S*V1).Q)
    
    @views if residual
        # The code for the residual case is adapted from 
        # DescriptorSystems.jl and written by Andreas Varga
        # https://github.com/andreasvarga/DescriptorSystems.jl/blob/dd144828c3615bea2d5b4977d7fc7f9677dfc9f8/src/order_reduction.jl#L622
        # with license https://github.com/andreasvarga/DescriptorSystems.jl/blob/main/LICENSE.md
        ONE = one(eltype(A))
        L = [Y Matrix(qr!(R'U2).Q)]
        Tr = [X Matrix(qr!(S*V2).Q)]
        amin = L'A*Tr
        bmin = L'B
        cmin = C*Tr
        Ar = amin[i1,i1]
        Er = L[:,i1]'*Tr[:,i1]
        Br = bmin[i1,:]
        Cr = cmin[:,i1]
        Dr = copy(D)
        LUF = lu!(amin[i2,i2])
        ldiv!(LUF,amin[i2,i1])
        ldiv!(LUF,bmin[i2,:])
        # apply state residualization formulas
        mul!(Dr,cmin[:,i2],bmin[i2,:],-ONE, ONE)
        mul!(Br,amin[i1,i2],bmin[i2,:],-ONE, ONE)
        mul!(Cr,cmin[:,i2],amin[i2,i1],-ONE, ONE)
        mul!(Ar,amin[i1,i2],amin[i2,i1],-ONE, ONE)
        # determine a standard reduced system
        SV = svd!(Er)
        di2 = Diagonal(1 ./sqrt.(SV.S))
        return ss(di2*SV.U'*Ar*SV.Vt'*di2, di2*(SV.U'*Br), (Cr*SV.Vt')*di2, Dr), Σ
    else
        L = (Y'X)\Y'
        T = X
        A = L*A*T
        B = L*B
        C = C*T
        D = D
    end
    ss(A,B,C,D, G.timeevol), Σ
end


"""
    controller_reduction_weight(P::ExtendedStateSpace, K)
Lemma 19.1 See Robust and Optimal Control Ch 19.1
"""
function controller_reduction_weight(P::ExtendedStateSpace, K)
    A, B1, B2, C1, C2, D11, D12, D21, D22 = ssdata_e(P)
    Ak,Bk,Ck,Dk = ssdata(K)
    R = lu(I - D22*Dk)
    R̃ = lu(I - Dk*D22)
    Aw = [
            A+B2*Dk*(R\C2) B2*(R̃\Ck)
            Bk*(R\C2)      Ak+Bk*D22*(R̃\Ck)
    ]
    Bw = [
        B2/R̃
        Bk*D22/R̃
    ]
    Cw = [R\C2   R\D22*Ck]
    Dw = D22/R̃
    ss(Aw, Bw, Cw, Dw)
end

"""
    controller_reduction(P::ExtendedStateSpace, K, r, out=true; kwargs...)

Minimize    ||(K-Kᵣ) W||∞ if out=false
            ||W (K-Kᵣ)||∞ if out=true
See Robust and Optimal Control Ch 19.1
out indicates if the weight will be applied as output or input weight.

This function expects a *positive feedback controller `K`.

This method corresponds to the methods labelled SW1/SW2(SPA) in 
Andreas Varga, "Controller Reduction Using Accuracy-Enhancing Methods"
SW1 is the default method, corresponding to `out=true`.

This method does not support unstable controllers. See the reference above for alternatives.
See also [`stab_unstab`](@ref) and [`baltrunc_unstab`](@ref).
"""
function controller_reduction(P::ExtendedStateSpace, K, r, out=true; kwargs...)
    W = controller_reduction_weight(P, K)
    if out
        frequency_weighted_reduction(K, W, 1, r; kwargs...)
    else
        frequency_weighted_reduction(K, 1, W, r; kwargs...)
    end
end

"""
    hsvd(sys::AbstractStateSpace)

Return the Hankel singular values of `sys`, computed as the eigenvalues of `QP`
Where `Q` and `P` are the Gramians of `sys`.
"""
function hsvd(sys::AbstractStateSpace)
    fun = ControlSystemsBase.isdiscrete(sys) ? MatrixEquations.plyapd : MatrixEquations.plyapc
    P = fun(sys.A, sys.B)
    Q = fun(sys.A', sys.C')
    e = svdvals(Q * P)
end

@deprecate minreal2 minreal

"""
    e = error_bound(hs)

Given a vector of Hankel singular values, return the theoretical error bound as a function of model order after balanced-truncation model reduction.
(twice sum of all the removed singular values).
`e[i]` gives you a bound of the error from removing modes `i:end`. For coprime controller reduction, you thus want to find the controller order
```
n = findlast(error_bound(hs) .> tolerated_ncfmargin)
```
"""
error_bound(hs) = [2reverse(cumsum(reverse(hs)))[1:end-1]; 0]

# slide 189 https://cscproxy.mpi-magdeburg.mpg.de/mpcsc/benner/talks/lecture-MOR.pdf
# This implementation works, but results in a complex-valued system.
# function model_reduction_irka(sys::AbstractStateSpace{Continuous}; n, tol = 1e-6)
#     A,B,C,D = ssdata(sys)
#     all(iszero, D) || error("Nonzero D not supported.")
#     sysr0, _ = baltrunc2(sys; n)
#     Ah, Bh, Ch = ssdata(sysr0)
#     e, T = eigen(Ah)
#     # _,T,E = blockdiagonalize(Ah)
#     # e = E.values
#     eold = 100e
#     k = 0
#     while maximum(abs.(e .- eold) ./ abs.(e)) > tol && k < 1000
#         Bt = T'Bh
#         Ct = Ch*T
#         Vv = map(1:n) do i
#             (-e[i]*I - A)\(B*Bt[i, :])
#         end
#         V = reduce(hcat, Vv)
#         Wv = map(1:n) do i
#             (-e[i]*I - A')\(C'*Ct[:,i])
#         end
#         W = reduce(hcat, Wv)
#         V = orthonormal(V) # not sure about those
#         W = orthonormal(W)
#         Ah = (W'V)\W'A*V
#         Bh = (W'V)\W'B
#         Ch = C*V
#         eold = e
#         k += 1
#         e, T = eigen(Ah)
#         # _,T,E = blockdiagonalize(Ah)
#         # e = E.values
#     end
#     sysr = ss(Ah, Bh, Ch, 0)
#     sysr, T
# end

# function orthonormal(X)
#     # U,S,V = svd(X, full=true)
#     # R = U*V'
#     Matrix(qr(X).Q) # |> real # this works sometimes
# end

"""
    controller_reduction_plot(G, K)

Plot the normalized-coprime margin ([`ncfmargin`](@ref)) as a function of controller order when [`baltrunc_coprime`](@ref) is used to reduce the controller.
Red, orange and green bands correspond to rules of thumb for bad, okay and good coprime uncertainty margins.
A value of 0 indicate an unstable closed loop.

If `G` is an ExtendedStateSpace system, a second plot will be shown indicating the \$H_∞\$ norm between inputs and performance outputs \$||T_{zw}||_\\infty\$ when the function [`controller_reduction`](@ref) is used to reduce the controller.

The order of the controller can safely be reduced as long as the normalized coprime margin remains sufficiently large. If the controller contains integrators, it may be advicable to protect the integrators from the reduction, e.g., if the controller is obtained using [`glover_mcfarlane`](@ref), perform the reduction on `info.Gs, info.Ks` rather than on `K`, and form `Kr` using the reduced `Ks`.

See [`glover_mcfarlane`](@ref) or [the docs](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/#Example-of-controller-reduction:) for an example.
"""
controller_reduction_plot
@userplot Controller_reduction_plot

@recipe function controller_reduction_plot(h::Controller_reduction_plot)
    G, K = h.args
    xguide --> "Model order"
    if G isa ExtendedStateSpace
        P = G
        G = system_mapping(G)
        norms = Float64[]
        failed = false
        for n = 1:K.nx-1
            try
                Kr, _ = controller_reduction(P, K, n, true)
                en = hinfnorm2(lft(P, Kr))[1]
                push!(norms, en)
            catch
                @info "Failed to reduce the controller to order using function controller_reduction"
                falied = true
            end
        end
        if !failed
            layout --> 2
            e = hinfnorm2(lft(P, K))[1]
            push!(norms, e)
            @series begin
                legend --> :topright
                yguide --> "||Tzw||∞"
                subplot --> 2
                label --> "\$||T_{zw}||_\\infty\$"
                linewidth --> 4
                norms
            end
        else
            layout --> 1
        end

    end

    subplot --> 1
    Kr, hs, info = baltrunc_coprime(K; n=1)
    margins = Float64[]
    for n = 1:K.nx-1
        Kr, _ = baltrunc_coprime(K, info; n)
        en = ncfmargin(G, Kr)[1]
        push!(margins, en)
    end
    e,_ = ncfmargin(G, K)
    push!(margins, e)
    yguide --> "Normalized coprime margin"
    legend --> :topleft
    @series begin
        label --> "Normalized coprime margin"
        linewidth --> 4
        margins
    end
    @series begin
        seriestype := :hline
        linestyle := :dash
        color --> :black
        label --> "ϵ = 0.25"
        fillbetween := (0.2, 0.3)
        fillcolor := :orange
        fillalpha := 0.12
        [0.25]
    end
    @series begin
        seriestype := :hline
        linestyle := :dash
        color --> :black
        label --> "ϵ = 0.25"
        fillbetween := (0.0, 0.2)
        fillcolor := :red
        fillalpha := 0.1
        primary := false
        [0.25]
    end
    @series begin
        seriestype := :hline
        linestyle := :dash
        color --> :black
        label --> "ϵ = 0.25"
        fillbetween := (0.3, 1)
        fillcolor := :green
        fillalpha := 0.12
        primary := false
        ylims := (0, 1)
        [0.25]
    end
end