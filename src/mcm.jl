function sys2Rn(fun, H::AbstractStateSpace)
    Ap,Bp,Cp,Dp = ssdata(H)
    PTF = ControlSystems.numeric_type(H)
    T,N,PT = MonteCarloMeasurements.particletypetuple(PTF)
    A = getindex.(Ap, 1)
    B = getindex.(Bp, 1)
    C = getindex.(Cp, 1)
    D = getindex.(Dp, 1)
    te = H.timeevol
    individuals = map(1:N) do i
        A .= getindex.(Ap, i)
        B .= getindex.(Bp, i)
        C .= getindex.(Cp, i)
        D .= getindex.(Dp, i)
        sys = ss(A,B,C,D,te)
        fun(sys)
    end
    MonteCarloMeasurements._finish_individuals(PT, Val(N), individuals, individuals[1])

end



function sys2sys(fun, H::AbstractStateSpace; allow_static=true, static_th = 10)
    Ap,Bp,Cp,Dp = ssdata(H)
    PTF = ControlSystems.numeric_type(H)
    ET,N,PT = MonteCarloMeasurements.particletypetuple(PTF)
    nx, nu, ny = H.nx, H.nu, H.ny
    if allow_static && size(Ap, 1) <= static_th
        individuals_s = static_kernel(fun, H, ET, PT, Val(nx), Val(nu), Val(ny), Val(nx*nx), Val(nx*nu), Val(nx*ny), Val(nu*ny))
        # individuals_s = static_kernel(fun, H, ET, PT, Val{nx}(), Val{nu}(), Val{ny}(), Val{nx*nx}(), Val{nx*nu}(), Val{nx*ny}(), Val{nu*ny}())
        return mergesystems(PTF, individuals_s)
    end
    A = getindex.(Ap, 1)
    B = getindex.(Bp, 1)
    C = getindex.(Cp, 1)
    D = getindex.(Dp, 1)
    te = H.timeevol
    individuals = map(1:N) do i
        A .= getindex.(Ap, i)
        B .= getindex.(Bp, i)
        C .= getindex.(Cp, i)
        D .= getindex.(Dp, i)
        sys = ss(A,B,C,D,te)
        fun(sys)
    end
    mergesystems(PTF, individuals)
end

function static_kernel(fun::F, H, ::Type{ET}, ::Type{PT}, ::Val{nx}, ::Val{nu}, ::Val{ny}, ::Val{nxnx}, ::Val{nxnu}, ::Val{nxny}, ::Val{nynu}) where {F,ET,PT,nx,nu,ny,nxnx,nxnu,nxny,nynu}
    Ap,Bp,Cp,Dp = ssdata(H)
    te = H.timeevol
    map(1:nparticles(Ap)) do pi
        A = SMatrix{nx, nx, ET, nxnx}(ntuple(i->Ap[i][pi], nxnx))
        B = SMatrix{nx, nu, ET, nxnu}(ntuple(i->Bp[i][pi], nxnu))
        C = SMatrix{ny, nx, ET, nxny}(ntuple(i->Cp[i][pi], nxny))
        D = SMatrix{ny, nu, ET, nynu}(ntuple(i->Dp[i][pi], nynu))
        sys = HeteroStateSpace(A,B,C,D,te)
        fun(sys)
    end
end

for PT in MonteCarloMeasurements.ParticleSymbols
    
    @eval function mergesystems(::Type{$PT{ET,N}},individuals::AbstractArray{<:StateSpace}) where {ET, N}
        PTF = $PT{ET, N} # full particle type
        s = individuals[1]
        nx,nu,ny = s.nx, s.nu, s.ny

        A::Matrix{PTF} = [
            PTF([individuals[pi].A[i,j] for pi in 1:N])
            for i in 1:nx, j in 1:nx
        ]
        B::Matrix{PTF} = [
            PTF([individuals[pi].B[i,j] for pi in 1:N])
            for i in 1:nx, j in 1:nu
        ]
        C::Matrix{PTF} = [
            PTF([individuals[pi].C[i,j] for pi in 1:N])
            for i in 1:ny, j in 1:nx
        ]
        D::Matrix{PTF} = [
            PTF([individuals[pi].D[i,j] for pi in 1:N])
            for i in 1:ny, j in 1:nu
        ]

        ss(A,B,C,D,s.timeevol)

    end
    @eval mergesystems(::Type{$PT{ET,N}},t::Vector{<:Tuple}) where {ET, N} = ntuple(i->mergesystems($PT{ET,N}, getindex.(t, i)), length(t[1]))
    @eval function mergesystems(::Type{$PT{ET,N}},a::Vector{<:AbstractArray}) where {ET, N}
        PTF = $PT{ET, N}
        nx = length(a[1])
        ret = [
            PTF([a[pi][i] for pi in 1:N])
            for i in 1:nx
        ]
        reshape(ret, size(a[1]))
    end
    # TODO: jag höll på med dessa, hade varit nice om tuple som returvärde hade funkat, försöker med att rekursivt kalla på mergesystems som då måste funka för arrays också
end



# for f in (:pole, :tzero)
#     @eval function ControlSystems.$f(s::TransferFunction{<:Any, <:ControlSystems.SisoRational{<:AbstractParticles}})
#         sys2Rn($f, s)
#     end
# end

for f in (:pole, :tzero)
    @eval function ControlSystems.$f(s::StateSpace{<:Any, <:AbstractParticles})
        sys2Rn($f, s)
    end
end

for f in (:_preprocess_for_freqresp, :minreal, :balreal)
    @eval function ControlSystems.$f(s::StateSpace{<:Any, <:AbstractParticles}, args...; kwargs...)
        sys2sys(s->$f(s, args...; kwargs...), s, allow_static=false)
    end
end

# function ControlSystems.balreal(s::StateSpace{<:Any, <:AbstractParticles}, args...; kwargs...)
#     sys2sys(s->balreal(s, args...; kwargs...)[1], s, allow_static=false)
# end

