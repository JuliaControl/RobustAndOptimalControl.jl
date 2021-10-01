using ControlSystems, RobustAndOptimalControl
plotly(size=(1850, 900), show=false)
##
w = 2pi*exp10.(LinRange(-2, 2, 400))
function flexible_servo_loop(; kp=1, kv=1, ki=1, Tf=0.1, Jm=1, Jl=10, c=1, k=10000)
    
    A = [
        0.0 1 0 0
        -k/Jm -c/Jm k/Jm c/Jm
        0 0 0 1
        k/Jl c/Jl -k/Jl -c/Jl
    ]
    B = [0, 1/Jm, 0, 0]
    C = I(4)
    D11 = D12 = D21 = D22 = 0
    nx,nu,ny = 4,1,1
    P = named_ss(ss(A,B,C,0), u=:u, x=[:qm, :qdm, :qa, :qda], y = [:qm, :qdm, :qa, :qda])
    
    K = let
        K = cascade_controller(; kp, kv, ki, Tf) # TODO: extend to take in several measurements instead of differentiating
        K = named_ss(K, u=:e, y = :Cu)
        # extended_controller(K)
    end 

    sumE = sumblock("e = r - y")

    u1 = [:u,  :y, :e]
    y1 = [:Cu, :qm, :e]
    w1 = [:r]
    G = connect([P, K, sumE]; y1, u1, w1)
    #F = ss([tf(0, [1, 0]); tf(1)])


    (; P, K, G)
end

function cascade_controller(; kp=1, kv=1, ki=1, Tf=0.1)
    ss(tf([kv, kp*kv + ki, kp*ki], [Tf, 1, 0]))
end

#
P,K,G = flexible_servo_loop(kp=0.1, kv=250, ki=15, Tf=0.001)

plot(
    # stepplot(c2d(G.sys, 0.001), 3),
    stepplot(G, 3),
    bodeplot(P, w, hz=true, plotphase=false),
    bodeplot(G, w, hz=true, plotphase=false),
    layout=(1,3),
    title="",
)