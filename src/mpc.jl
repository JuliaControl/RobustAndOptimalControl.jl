using ControlSystemsBase, Dates, LinearAlgebra, Plots, StaticArrays
const V32 = Vector{Float32}
const M32 = Matrix{Float32}
const V64 = Vector{Float64}
const M64 = Matrix{Float64}
const A32 = Union{Array{Float32}, <:SubArray{Float32, 1, <:Array{Float32}}, <:SubArray{Float32, 2, <:Array{Float32}}}

cd(@__DIR__)

function compile_lib(dir::String)
    print("Compiling library to ", dir, "\n")
    build_dir = joinpath(dir, "build")
    if !isdir(build_dir)
        mkdir(build_dir)
    end
    run(`cmake -S $dir -B $build_dir`)
    run(`cmake --build $build_dir`)
    return true
end

get_os_extension() = Sys.isapple() ? ".dylib" : Sys.iswindows() ? ".dll" : ".so"


tinympc_dir = "/home/fredrikb/repos/tinympc-julia/tinympc/TinyMPC"  # Path to the TinyMPC directory (C code)
compile_lib(tinympc_dir)  # Compile the C code into a shared library

tinympc = joinpath(tinympc_dir, "build","src","tinympc","libtinympcShared")*get_os_extension()  # Path to the compiled library
@assert isfile(tinympc)  # Check that the library exists

struct MPCController2
    x::M32
    u::M32
    library_path::String
    library::Ptr{Cvoid}
    set_x0_ptr::Ptr{Cvoid}
    set_xref_ptr::Ptr{Cvoid}
    call_tiny_solve_ptr::Ptr{Cvoid}
    get_u_ptr::Ptr{Cvoid}
    get_x_ptr::Ptr{Cvoid}
    t::Vector{Float64}
end

Base.size(c::MPCController2) = size(c.u), size(c.x)

function MPCController2(sys, Q1::M64, Q2::M64, N::Integer;
    x_min::V64 = fill(-1e6, sys.nx*N),  # state constraints
    x_max::V64 = fill(1e6, sys.nx*N),  # state constraints
    u_min::V64 = fill(-1e6, sys.nu*(N-1)),  # input constraints
    u_max::V64 = fill(1e6, sys.nu*(N-1)),  # input constraints
    rho::Float64 = 0.1,
    abs_pri_tol::Float64 = 1.0e-3,  # absolute primal tolerance
    abs_dual_tol::Float64 = 1.0e-3,  # absolute dual tolerance
    max_iter::Integer = 10000,  # maximum number of iterations
    check_termination::Integer = 2,  # whether to check termination and period
    output_dir = "generated_code_$(now())",  # Path to the generated code
    verbose::Integer = true,
    compile = true,
)

    (; A, B, nx, nu) = sys
    isdiag(Q1) || throw(ArgumentError("Q1 must be diagonal"))
    isdiag(Q2) || throw(ArgumentError("Q2 must be diagonal"))
    size(Q1) == (nx, nx) || throw(ArgumentError("Q1 must have size nx x nx"))
    size(Q2) == (nu, nu) || throw(ArgumentError("Q2 must have size nu x nu"))
    if length(x_min) == nx
        x_min = repeat(x_min, N)
    end
    if length(x_max) == nx
        x_max = repeat(x_max, N)
    end
    if length(u_min) == nu
        u_min = repeat(u_min, N-1)
    end
    if length(u_max) == nu
        u_max = repeat(u_max, N-1)
    end
    length(x_min) == sys.nx*N || throw(ArgumentError("x_min must have length nx*N"))
    length(x_max) == sys.nx*N || throw(ArgumentError("x_max must have length nx*N"))
    length(u_min) == sys.nu*(N-1) || throw(ArgumentError("u_min must have length nu*(N-1)"))
    length(u_max) == sys.nu*(N-1) || throw(ArgumentError("u_max must have length nu*(N-1)"))

    if !(eltype(A) <: Float64)
        A, B = convert(Matrix{Float64}, A), convert(Matrix{Float64}, B)
    end
    @ccall tinympc.tiny_codegen(Cint(nx)::Cint, Cint(nu)::Cint, Cint(N)::Cint, A::Ptr{Float64}, B::Ptr{Float64}, diag(Q1)::Ptr{Float64}, diag(Q2)::Ptr{Float64}, x_min::Ptr{Float64}, x_max::Ptr{Float64}, u_min::Ptr{Float64}, u_max::Ptr{Float64}, rho::Float64, abs_pri_tol::Float64, abs_dual_tol::Float64, Cint(max_iter)::Cint, Cint(check_termination)::Cint, Cint(verbose)::Cint, tinympc_dir::Ptr{UInt8}, output_dir::Ptr{UInt8})::Cint

    library_path = joinpath(output_dir,"build","tinympc","libtinympcShared")*get_os_extension()
    if compile
        compile_lib(output_dir)
        library = Libc.dlopen(library_path)
        set_x0_ptr = Libc.dlsym(library, :set_x0)
        set_xref_ptr = Libc.dlsym(library, :set_xref)
        call_tiny_solve_ptr = Libc.dlsym(library, :call_tiny_solve)
        get_u_ptr = Libc.dlsym(library, :get_u)
        get_x_ptr = Libc.dlsym(library, :get_x)
        ccall(set_xref_ptr, Cvoid, (Ptr{Float32}, Cint), zeros(Float32, N*nx), Cint(0))
    else
        library = Ptr{Cvoid}(0)
        set_x0_ptr = Ptr{Cvoid}(0)
        set_xref_ptr = Ptr{Cvoid}(0)
        call_tiny_solve_ptr = Ptr{Cvoid}(0)
        get_u_ptr = Ptr{Cvoid}(0)
        get_x_ptr = Ptr{Cvoid}(0)
    end

    MPCController2(
        zeros(Float32, sys.nx, N),
        zeros(Float32, sys.nu, (N-1)),
        library_path,
        library,
        set_x0_ptr,
        set_xref_ptr,
        call_tiny_solve_ptr,
        get_u_ptr,
        get_x_ptr,
        Float64[],
    )

end

function (controller::MPCController2)(x0::A32, r::Union{Nothing, A32}=nothing; verbose=false)
    (; set_x0_ptr,
        set_xref_ptr,
        call_tiny_solve_ptr,
        get_u_ptr,
        get_x_ptr) = controller

    if r !== nothing
        if length(r) == size(controller.x, 1)
            r = repeat(r, size(controller.x, 2))
        end
        length(r) == length(controller.x) || throw(ArgumentError("r must have the same length as the state vector"))
        ccall(set_xref_ptr, Cvoid, (Ptr{Float32}, Cint), r, Cint(0))
    end

    ccall(set_x0_ptr, Cvoid, (Ptr{Float32}, Cint), x0, Cint(0))

    t = @elapsed ccall(call_tiny_solve_ptr, Cvoid, (Cint,), Cint(verbose))
    push!(controller.t, t)
    ccall(get_u_ptr, Cvoid, (Ptr{Float32}, Cint), controller.u, 0)
    ccall(get_x_ptr, Cvoid, (Ptr{Float32}, Cint), controller.x, 0)

    (; controller.u, controller.x, t)
end

##
using RobustAndOptimalControl
Ts = 0.05
N = 100
Pc = DemoSystems.double_mass_model(c0=0.01,c1=0.01,c2=0.01, outputs=1)
P = c2d(Pc, Ts)
P = add_input_integrator(P, 1; ϵ=1e-3)
Q1 = Matrix(1.0*I(P.nx))
Q2 = Matrix(0.1*I(P.nu))


u_min = fill(-250.0, P.nu)
u_max = fill(250.0, P.nu)
x_max = [1000.0, 1000, 1000, 1000, 10]
x_min = -x_max

##

Lmpc = MPCController2(P, Q1, Q2, N;
    x_min,
    x_max,
    u_min,
    u_max,
    rho = 1.0,
    max_iter = 3500,
    abs_pri_tol = 1.0e-3,
    abs_dual_tol = 1.0e-3,
    check_termination = 2,

)
L = lqr(P, Q1, Q2)

# x0 = zeros(Float32, P.nx)
x0 = Float32[20, 0, 0, 0, 0]
r = reduce(hcat, fill(zeros(Float32, P.nx), N))

# u, x, t = Lmpc(x0, r, verbose=false)
# plot((P.C*x)', layout=(2,1)); plot!(u', sp=2)

lqr_fun = (x,t)->clamp.(-L*x, u_min[1], u_max[1])
mpc_closure =function (P, nu::Val{NU} = Val(P.nu)) where NU
    x_F32 = zeros(Float32, P.nx)
    u_F32 = zeros(Float32, P.nu)
    let Lmpc = Lmpc, u_min = u_min, u_max = u_max
        function (x,t)
            x_F32 .= x
            res = Lmpc(x_F32; verbose=false)
            @views u_F32 .= clamp.(
                res.u[1:NU, 1],
                u_min[1:NU, 1],
                u_max[1:NU, 1]
            )
            SVector{NU}(u_F32)
        end
    end

end

mpc_fun = mpc_closure(P)

@time "lqr" res_lqr = lsim(P, lqr_fun, 5; x0);
@time "mpc" res_mpc = lsim(P, mpc_fun, 5; x0);

cost_lqr = dot(res_lqr.x, Q1, res_lqr.x) + dot(res_lqr.u, Q2, res_lqr.u)
cost_mpc = dot(res_mpc.x, Q1, res_mpc.x) + dot(res_mpc.u, Q2, res_mpc.u)

maximum_constraint_violation = max(
    maximum(res_mpc.x .- x_max),
    maximum(x_min .- res_mpc.x),
    maximum(res_mpc.u .- u_max),
    maximum(u_min .- res_mpc.u)
)


fig1 = plot(res_lqr, label="LQR $cost_lqr", plotu=true, plotx=false, size=(800, 1200), margin=5Plots.mm)
plot!(res_mpc, label="MPC $cost_mpc", plotu=true, plotx=false)
fig2 = scatter(1e3 .* Lmpc.t, title="TimyMPC execution time [ms]", label=false)
plot(fig1, fig2, layout=(1,2), size=(1200, 1200))
display(current())

