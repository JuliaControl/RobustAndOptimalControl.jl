#=
This example builds a control system for a servo motor that is connected to a load through a flexible transmission. The servo is controlled using a cascade controller with a PI controller for the velocity and a P controller for the position.
=#
using ControlSystemsBase, RobustAndOptimalControl, Plots

function flexible_servo_loop(;
    kp=1, # proportional gain for position loop
    kv=1, # proportional gain for velocity loop (derivative gain)
    ki=1, # integral gain for the velocity loop
    Tf=0.1, # filter time constant
    Jm=1,   # inertia of motor
    Jl=10,  # inertia of load
    c=1,    # damping
    k=10000 # spring stiffness
)
    
    A = [
        0.0 1 0 0
        -k/Jm -c/Jm k/Jm c/Jm
        0 0 0 1
        k/Jl c/Jl -k/Jl -c/Jl
    ]
    B = [0, 1/Jm, 0, 0]
    C = I(4) # All states are measurable (not realistic, only included to show Bode plots)
    D11 = D12 = D21 = D22 = 0
    nx,nu,ny = 4,1,1
    P = named_ss(ss(A,B,C,0), u=:u, x=[:qm, :qdm, :qa, :qda], y = [:qm, :qdm, :qa, :qda])
    
    K = let
        K = cascade_controller(; kp, kv, ki, Tf) # TODO: extend to take in several measurements instead of differentiating
        K = named_ss(K, u=:e, y = :Cu)
        # extended_controller(K)
    end 

    sumE = sumblock("e = r - y") # Form the summing node that computes the control error


    connections = [
        :Cu => :u
        :qm => :y
        :e => :e
    ]
    w1 = [:r] # r is an external input
    G = connect([P, K, sumE], connections; w1)

    # Alternative calling convention:
    # u1 = [:u,  :y, :e]   # inputs
    # y1 = [:Cu, :qm, :e]  # outputs
    # G = connect([P, K, sumE]; u1, y1, w1)

    (; P, K, G)
end

function cascade_controller(; kp=1, kv=1, ki=1, Tf=0.1)
    ss(tf([kv, kp*kv + ki, kp*ki], [Tf, 1, 0]))
end

#
P,K,G = flexible_servo_loop(kp=0.1, kv=250, ki=15, Tf=0.001)
w = 2pi*exp10.(LinRange(-2, 2, 400)) # Frequency vector

plot(
    # plot(step(c2d(G, 0.001), 3)),
    plot(step(G, 3, method=:zoh)),
    bodeplot(P, w, hz=true, plotphase=false, title="P"),
    bodeplot(G, w, hz=true, plotphase=false, title="Closed loop", legend=:bottomleft),
    layout = (1,3),
    title  = "",
    size   = (1000, 1000),
)
