using RobustAndOptimalControl, ControlSystems, Plots
P = named_ss(ss(tf(1, [1, 0.1, 1])), :P)
We = named_ss(makeweight(10, 1, 0.1),  :We, u=:y, y=:e)  # We care a lot about the error on low frequencies
Wu = named_ss(0.01makeweight(1e-3, 10, 10), :Wu, u=:Wu, y=:uw) # Above ω=100, we want to limit the control effort
# Wd = named_ss(makeweight(1, 1, 1e-3),  :Wd, u=:do, y=:d) # d is low frequency
Wd = named_ss(ss(1), :Wd, u=:do, y=:d)

sumP = sumblock("y = Py + d")


split_u = splitter(:u, 2)


connections = [
    :u1 => :Wu # splitter to input of Wu
    :u2 => :Pu # splitter to input of P
    :Py => :Py # P output to first input of sumblock
    :d => :d   # output of Wd to second input of sumblock
    :y => :y   # output of sumblock to input of We
];

w1 = [ # External inputs
    :do, :u
]
z1 = [ # External outputs
    :e, :uw, :y
];

G = connect([P,We,Wu,Wd,sumP,split_u], connections; z1, w1)

Gsyn = partition(G, u = [:u], y = [:y]) # You can provide either u or w, and either y or z
K, γ, info = hinfsynthesize(Gsyn, γrel=1.001, interval = (0.1, 20), transform=false)


Gsyn2 = hinfpartition(P, We.sys, Wu.sys, [])
K2, γ2 = hinfsynthesize(Gsyn2, γrel=1.001, interval = (0.1, 20), transform=false)

@test γ ≈ 0.3148 atol=1e-2 # value by slicot
@test γ ≈ γ2 atol=1e-3

@test hinfnorm2(lft(Gsyn, K))[1] ≈ γ atol=1e-2
@test hinfnorm2(lft(Gsyn2, K2))[1] ≈ γ2 atol=1e-2


@test hinfnorm2(K+K2)[1] < 1e-6


Pcl, S, CS, T = hinfsignals(Gsyn, P, K)

@test hinfnorm2(Pcl)[1] ≈ γ atol=1e-2
@test hinfnorm2(lft(Gsyn, K))[1] ≈ γ atol=1e-2
@test hinfnorm2(lft(Gsyn2, K2))[1] ≈ γ2 atol=1e-2


# Controller by slicot, test if hinfnorm is the same
Km = let
    tempA = [0.0 1.0 0.0 0.0; -227.5648 -32.5972 -703.6977 -117.9357; 0.0 0.0 -0.1 0.0; -226.5648 -32.4972 -703.6977 -217.4345]
    tempB = [-0.0; -0.0; 1.0; 0.0;;]
    tempC = [-226.5648 -32.4972 -703.6977 -117.9357]
    tempD = [0.0;;]
    ss(tempA, tempB, tempC, tempD)
end

@test hinfnorm2(lft(Gsyn, Km))[1] ≈ 0.3148 atol=1e-3
@test hinfnorm(lft(Gsyn, Km))[1] ≈ 0.3148 atol=1e-3



