show_construction(sys::LTISystem; kwargs...) = show_construction(stdout, sys; kwargs...)

"""
    show_construction([io::IO,] sys::LTISystem; name = "temp", letb = true)

Print code to `io` that reconstructs `sys`.
- `letb`: If true, the code is surrounded by a let block.

```@example
julia> sys = ss(tf(1, [1, 1]))
StateSpace{Continuous, Float64}
A = 
 -1.0
B = 
 1.0
C = 
 1.0
D = 
 0.0

Continuous-time state-space model

julia> show_construction(sys, name="Jörgen")
Jörgen = let
    JörgenA = [-1.0;;]
    JörgenB = [1.0;;]
    JörgenC = [1.0;;]
    JörgenD = [0.0;;]
    ss(JörgenA, JörgenB, JörgenC, JörgenD)
end
```
"""
function show_construction(io::IO, sys::LTISystem; name = "temp", letb = true)
    # sys = StateSpace(sys)
    letb && println(io, "$name = let")
    prestr = letb ? "    " : "" 
    println(io, prestr*"$(name)A = ", sys.A)
    println(io, prestr*"$(name)B = ", sys.B)
    println(io, prestr*"$(name)C = ", sys.C)
    println(io, prestr*"$(name)D = ", sys.D)
    letb || print(io, "$name = ")
    if isdiscrete(sys)
        println(io, prestr*"ss($(name)A, $(name)B, $(name)C, $(name)D, $(sys.Ts))")
    else
        println(io, prestr*"ss($(name)A, $(name)B, $(name)C, $(name)D)")
    end
    letb && println(io, "end")
    nothing
end

function show_construction(io::IO, sys::NamedStateSpace; name = "temp", letb = true)
    # sys = StateSpace(sys)
    letb && println(io, "$name = let")
    prestr = letb ? "    " : "" 
    println(io, prestr*"$(name)A = ", sys.A)
    println(io, prestr*"$(name)B = ", sys.B)
    println(io, prestr*"$(name)C = ", sys.C)
    println(io, prestr*"$(name)D = ", sys.D)
    letb || print(io, "$name = ")
    if isdiscrete(sys)
        println(io, prestr*"named_ss(ss($(name)A, $(name)B, $(name)C, $(name)D), $(sys.Ts), x=$(sys.x), u=$(sys.u), y=$(sys.y))")
    else
        println(io, prestr*"named_ss(ss($(name)A, $(name)B, $(name)C, $(name)D), x=$(sys.x), u=$(sys.u), y=$(sys.y))")
    end
    letb && println(io, "end")
    nothing
end

function Base.vec(sys::LTISystem)
    [vec(sys.A); vec(sys.B); vec(sys.C); vec(sys.D)]
end


"""
    vec2sys(v::AbstractArray, ny::Int, nu::Int, ts = nothing)

Create a statespace system from the parameters
```julia
v = vec(sys) = [vec(sys.A); vec(sys.B); vec(sys.C); vec(sys.D)]
```
Use `vec(sys)` to create `v`.

This can be useful in order to convert to and from vectors for, e.g., optimization.

```@example
julia> sys  = ss(tf(1, [1, 1]))
StateSpace{Continuous, Float64}
A = 
 -1.0
B = 
 1.0
C = 
 1.0
D = 
 0.0

Continuous-time state-space model

julia> v    = vec(sys)
4-element Vector{Float64}:
 -1.0
  1.0
  1.0
  0.0

julia> sys2 = vec2sys(v, sys.ny, sys.nu)
StateSpace{Continuous, Float64}
A = 
 -1.0
B = 
 1.0
C = 
 1.0
D = 
 0.0

Continuous-time state-space model
```
"""
function vec2sys(v::AbstractArray, ny::Int, nu::Int, ts=nothing)
    n = length(v)
    p = (ny+nu)
    nx = Int(-p/2 + sqrt(p^2 - 4nu*ny + 4n)/2)
    @assert n == nx^2 + nx*nu + ny*nx + ny*nu
    ai = (1:nx^2)
    bi = (1:nx*nu) .+ ai[end]
    ci = (1:nx*ny) .+ bi[end]
    di = (1:nu*ny) .+ ci[end]
    A = reshape(v[ai], nx, nx)
    B = reshape(v[bi], nx, nu)
    C = reshape(v[ci], ny, nx)
    D = reshape(v[di], ny, nu)
    ts === nothing ? ss(A, B, C, D) : ss(A, B, C, D, ts)
end

# sys2vec = @(sys) [
#         size(sys.A,1);
#         size(sys.B,2);
#         size(sys.C,1);
#         sys.A(:);
#         sys.B(:);
#         sys.C(:);
#         sys.D(:)
#     ]
# function vec2sys(v::AbstractArray, ts = nothing; kwargs...)
#     nx = Int(v[1])
#     nu = Int(v[2])
#     ny = Int(v[3])
#     ai = (1:nx^2) .+ 3
#     bi = (1:nx*nu) .+ ai[end]
#     ci = (1:nx*ny) .+ bi[end]
#     di = (1:nu*ny) .+ ci[end]
#     A = reshape(v[ai], nx, nx)
#     B = reshape(v[bi], nx, nu)
#     C = reshape(v[ci], ny, nx)
#     D = reshape(v[di], ny, nu)
#     sys = ts === nothing ? ss(A, B, C, D) : ss(A, B, C, D, ts)
#     show_construction(sys; kwargs...)
#     sys
# end