using ControlSystems, RobustAndOptimalControl, MonteCarloMeasurements
##

a = 10
P = [
        tf([1,-a^2], [1, 0, a^2]) tf([a, a], [1, 0, a^2])
        -tf([a, a], [1, 0, a^2]) tf([1,-a^2], [1, 0, a^2])
    ]
P = minreal(ss(P))
K = ss(1.0I(2))


ny,nu = size(P)
sys2 = ss(I(ny+nu)) # this formulation makes sense if sys2 is I + a*δ 

# sys1 = ControlSystems.append(P,K)
# z1 = 1:ny
# w1 = (1:ny) .+ ny
# M = feedback(sys1, sys2; Z1 = z1, W1 = w1)

# M = feedback(P, K, W2 = :, Z2 = :) # output everything
# feedback(M, ss(0*I(ny+nu)))

# 1. form system that exposes all inputs and outputs but also has feedback
M = feedback(P, K, W2 = :, Z2 = :)
@test poles(M) ≈ poles(feedback(P*K))
@test size(M) == (ny+nu, ny+nu)
# @test minreal(feedback(M, 0*sys2)) ≈ M


w = 2π .* exp10.(LinRange(-2, 2, 300))
# @time bisect_a(P, K, w)
##

# break at input (pass outputs through)
a = bisect_a(P, K, w; Z = [], tol=1e-4)
@test minimum(a) < 0.0998

# break at output (pass inputs through)
a = bisect_a(P, K, w; W = [], tol=1e-4)
@test minimum(a) < 0.0998

# break at both input and output
a = bisect_a(P, K, w; tol=1e-4)
au = bisect_a(P, K, w; tol=1e-4, N=640, upper=true)
@test minimum(a) < 0.0499


plot(w, a, xscale=:log10, xlabel="Frequency", ylims=(0,3))
plot!(w, au, xscale=:log10, xlabel="Frequency", ylims=(0,3))


##
w = 2π .* exp10.(LinRange(-2, 2, 300))

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


a = bisect_a(L3, ss(I(3), L3.Ts), w; tol=2e-3)
au = bisect_a(L3, ss(I(3), L3.Ts), w; tol=2e-3, upper=true, N=256)
plot(w, a, xscale=:log10, xlabel="Frequency", ylims=(0,Inf))
plot!(w, au, xscale=:log10, xlabel="Frequency", ylims=(0,Inf))


dm = diskmargin(L3, 0, w)
plot(w, dm[1,:])


dm = diskmargin(L3, 0, 4.05)
@test dm[1].α ≈ 0.794418036911981 rtol=1e-3


dm = diskmargin(L3, ss(I(3), L3.Ts), 0, w)


## Test diskmargin with particles in the system

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]*(1 + 0.1*Particles(32))
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


unsafe_comparisons(true)
dm = diskmargin(L3, 0, 4.05)




L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


P = L3
K = ss(I(3), L3.timeevol)
w = 2π .* exp10.(LinRange(-2, 2, 300))
# break at input (pass outputs through)
M,_ = get_M(P, K, w; Z = :)
M0 = M[:,:,100]



@time mu = [mussv(M) for i = 1:10]
# mum = minimum(mu)
plot(w, mu, xscale=:log10)

a = bisect_a(P, K, w; Z = [], tol=1e-4)
au = bisect_a(P, K, w; Z = [], tol=1e-4, upper=true, N=2560)
plot!(w, inv.(a), xscale=:log10)
plot!(w, inv.(au), xscale=:log10)

## ==================================================


C = K
dm = diskmargin(L3, 1.0*K, 0, w;)
plot(dm.simultaneous) # TODO: verkar ge fel svar, ungefär samma form som matlab men inte rätt värden

plot(dm.input)





## 
# NOTE: SISO och loop at a time blir rätt, men inte simultaneous. mussv verkar rätt, så kan vara fel på get_M
a = [-0.2 10;-10 -0.2]; b = I(2); c = [1 8;-10 1];
P = ss(a,b,c,0);
K = ss([1 -2;0 1]);
dm = diskmargin(K*P) # disk margins at plant inputs
dm = diskmargin(P*K); # disk margins at plant outputs
MMIO = diskmargin(P,K,w)

plot(MMIO.simultaneous)

w = 2π .* exp10.(LinRange(-3, 3, 300))

##
MMO = diskmargin(P,K,0,w)

plot(MMO.simultaneous, lab="simultaneous")
plot!(MMO.simultaneous_output, lab="output") # denna är rätt för små frekvenser men 2x fel för höga
plot!(MMO.simultaneous_input, lab="input") # samma för denna




L = K*P
M0 = permutedims(freqresp(feedback(L), w), (2,3,1))
mu = mussv(M0)
imu = inv.(mussv(M0))
simultaneous = [Diskmargin(imu; ω0 = w, L) for (imu, w) in zip(imu,w)]

plot(w, mu)


##
a0 = 10 *(1 + δ(8))
a = [-0.2 a0;-a0 -0.2]; b = I(2); c = [1 8;-a0 1];
P = ss(a,b,c,0);
K = ss([1 -2;0 1]);

w = 2π .* exp10.(LinRange(-1, log10(500), 500))
freqresp(P, w)

sys = P
s = 0.1im

R\sys.B
