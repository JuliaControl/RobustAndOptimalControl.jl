using ControlSystems, RobustAndOptimalControl
function flexible_servo_loop(; kp=1, kv=1, ki=1, Tf=0.1, Jm=1, Jl=10, c=1, k=10000)
    
    A = [
        0.0 1 0 0
        -k/Jm -c/Jm k/Jm c/Jm
        0 0 0 1
        k/Jl c/Jl -k/Jl -c/Jl
    ]
    B = [0, 1/Jm, 0, 0]
    C2 = [
        1 0 0 0
        # 0 1 0 0
    ]
    C1 = I(4)
    D11 = D12 = D21 = D22 = 0
    nx,nu,ny = 4,1,1
    P = ExtendedStateSpace(A,zeros(nx,0),B,C1,C2)
    
    K = let
        K = cascade_controller(; kp, kv, ki, Tf) # TODO: extend to take in several measurements instead of differentiating
        nx,nu,ny = K.nx, K.nu, K.ny
        A,B,C,D = ssdata(K)
        ss(A, B, -B, zeros(0,nx), C, D21=D, D22=-D)
    end

    G = lft(P,K)
    #F = ss([tf(0, [1, 0]); tf(1)])

    (; P, K, G)
end

function cascade_controller(; kp=1, kv=1, ki=1, Tf=0.1)
    ss(tf([kv, kp*kv + ki, kp*ki], [Tf, 1, 0]))
end

##
P,K,G = flexible_servo_loop(kp=0.1, kv=250, ki=150, Tf=0.001)

plot(
    stepplot(ss(G), 3),
    bodeplot(ss(P), 2pi*exp10.(LinRange(-2, 2, 400)), hz=true, plotphase=false),
    bodeplot(ss(G), 2pi*exp10.(LinRange(-2, 2, 400)), hz=true, plotphase=false),
    layout=(1,3),
    title="",
)