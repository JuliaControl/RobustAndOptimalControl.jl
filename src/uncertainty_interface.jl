
function δp(a = 1, N=32)
    r = Particles(N, Uniform(0, a))
    θ = Particles(N, Uniform(0, 2pi))
    r*cis(θ)
end

function LinearAlgebra.det(D::AbstractMatrix{<:Complex{<:AbstractParticles}})
    D0 = similar(D, ComplexF64)
    parts = map(1:nparticles(D[1].re)) do i
        for j in eachindex(D0)
            D0[j] = Complex(D[j].re.particles[i], D[j].im.particles[i])
        end
        det(D0)
    end
    Complex(StaticParticles(getfield.(parts, :re)), StaticParticles(getfield.(parts, :im)))
end

function ControlSystems.tzeros(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractMatrix{T}, D::AbstractMatrix{T}) where T <: AbstractParticles
    bymap(tzeros, A, B, C, D)
end


using MonteCarloMeasurements: vecindex
function sys_from_particles(P, i)
    A,B,C,D = ssdata(P)
    ss(vecindex(A, i), vecindex(B, i), vecindex(C, i), vecindex(D, i))
end


function any0det(D::Matrix{<:Complex{<:AbstractParticles}})
    D0 = similar(D, ComplexF64)
    maxre = maxim = -1
    minre = minim = 1
    for i = 1:nparticles(D[1].re)
        for j in eachindex(D0)
            D0[j] = Complex(D[j].re.particles[i], D[j].im.particles[i])
        end
        d = det(D0)
        maxre = max(maxre, d.re)
        minre = min(minre, d.re)
        maxim = max(maxim, d.im)
        minim = min(minim, d.im)
        if maxre > 0 && minre < 0 && maxim > 0 && minim < 0
            return true
        end
    end
    false
end
