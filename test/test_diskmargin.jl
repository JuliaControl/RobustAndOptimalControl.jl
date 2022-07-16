using ControlSystems, RobustAndOptimalControl, MonteCarloMeasurements
# using RobustAndOptimalControl: bisect_a

# Example from the diskmargin paper

L = tf(25, [1,10,10,10])
dm = diskmargin(L, 0)
@test dm.ω0 ≈ 1.94   atol=0.02
@test dm.γmin ≈ 0.63    atol=0.02
@test dm.γmax ≈ 1.59    atol=0.02
@test dm.α ≈ 0.46   atol=0.02
show(dm)
plot(dm)
nyquistplot(L)
plot!(dm, nyquist=true)
plot!(Disk(dm), nyquist=true)

@test dm.alpha == dm.α
@test length(dm.gainmargin) == 2
@test dm.phasemargin == dm.ϕm
@test :gainmargin ∈ propertynames(dm)



## Frequency-dependent margin
w = exp10.(LinRange(-2, 2, 500))
dms = diskmargin(L, 0, w)
plot(dms)

##
s = tf("s")
L = 6.25*(s + 3)*(s + 5) / (s*(s + 1)^2 *(s^2 + 0.18s + 100))

## αmax > 2
dm = diskmargin(L, 0, 200)
@test dm.γmax < dm.γmin
@test dm.γmin ≈ 0 atol = 1e-6
@test dm.ϕm ≈ 90 atol = 1e-4


w = exp10.(LinRange(-1, 2, 500))
dms = diskmargin(L, 0, w)
plot(dms) 



# using IntervalArithmetic
# δ(a=1) = Complex(-a..a, -a..a)
# Δ(n, a) = diagm([δ(a) for _ in 1:n])
# M = [0 1; -0.1 -0.1]
# D = Δ(2, 1)
# 0 ∈ det(I-M*D)


## Loop at a time
a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)
K = ss(I(2))

Li = K*P
Lo = P*K

@test tf(minreal(broken_feedback(Li, 1))) ≈ tf(1, [1, 0])
@test tf(minreal(broken_feedback(Lo, 1))) ≈ tf(1, [1, 0])

dm = loop_diskmargin(Li)[1]
@test dm.α ≈ 2
@test dm.γmax > 1e10
@test dm.ϕm ≈ 90

dmm = loop_diskmargin(P, K)
dmm.input[1].α == dm.α


##
L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end

w = 2π .* exp10.(LinRange(-2, 2, 300))

dm = loop_diskmargin(L3, 0, 4.05)
@test dm[1].α ≈ 0.794418036911981 rtol=1e-3

dm = diskmargin(L3, ss(I(3), L3.Ts), 0, w)
plot(dm.simultaneous)
plot!(dm.simultaneous_input)
plot!(dm.simultaneous_output)

plot(dm.input)
plot!(dm.output)

## MIMO

a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)
K = ss(I(2))

w = 2π .* exp10.(LinRange(-2, 2, 300))
##

# break at input (pass outputs through)
# a = bisect_a(P, K, w; Z = [], tol=1e-4)

@test diskmargin(K*P).α ≈ 0.0998 atol=0.002 # from DM paper

# @test minimum(a) < 0.0998 

# break at output (pass inputs through)
# a = bisect_a(P, K, w; W = [], tol=1e-4)
# @test minimum(a) < 0.0998

# break at both input and output
# a = bisect_a(P, K, w; tol=1e-4)
# au = bisect_a(P, K, w; tol=1e-4, N=640, upper=true)
@test sim_diskmargin(P, K).α ≈ 0.0499 atol=0.002

# ar = bisect_a(P, K, w; tol=1e-4, δ=δr)
# aur = bisect_a(P, K, w; tol=1e-4, N=640, upper=true, δ=δr)

# dm = diskmargin(P, K, 0, w)
# adm = [dm.α for dm in dm.simultaneous]


# if isinteractive()
#     plot(w,  a,   xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Lower Complex")
#     plot!(w, au,  xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Upper Complex")
#     plot!(w, ar,  xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Lower Real")
#     plot!(w, aur, xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Upper Real")
#     plot!(w, adm, xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="μ")
#     display(current())
# end


##
w = 2π .* exp10.(LinRange(-2, 2, 300))

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end

dm = diskmargin(L3, 0, w)
plot(dm)

## Test diskmargin with particles in the system

L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]*(1 + 0.1*Particles(32))
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


unsafe_comparisons(true)
dm = loop_diskmargin(L3, 0, 4.05)
@test dm[1].α isa Particles
@test pmean(dm[1].α) ≈ 0.794418036911981 rtol=0.01


L3 = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    tempD = [0.0 0.0 0.0; 34875.36444283988 0.0 0.0; 48304.01940122544 0.0 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end


P = L3
K = ss(1.0I(3), L3.timeevol)

w = [0.01, 0.1, 1, 10, 100]
mu = 1e5*[2.961061403927230
2.961017564832289
0.029884614923336
0.029880921297303
0.000453111209111
0.000452684368322
0.000033811562301
0.000033712186847
0.000003489834242
0.000003488286208]

mu2 = structured_singular_value(freqresp(L3, w))
@test mu2 ≈ mu[2:2:end] rtol=1e-6


##
w = exp10.(LinRange(-2, 2, 300))
K = ss(I(3), L3.timeevol)
L = feedback(L3, K)
mu = structured_singular_value(freqresp(L, w))
plot(mu)
@test mu[1] ≈ 1 rtol=1e-3
@test maximum(mu) ≈ 2.503148529400597 rtol=1e-3
@test mu[end] ≈ 0.397718223553195 rtol=0.03

##
w = 2π .* exp10.(LinRange(-2, 2, 300))
# break at input (pass outputs through)
# M,D = RobustAndOptimalControl.get_M(P, K, w; W = [])
# M,D = RobustAndOptimalControl.get_M(L3, w)


# Test stability of the computation
# @time mu = [structured_singular_value(M) for i = 1:10]
# plot(w, mu, xscale=:log10)


## compare structured_singular_value and bisect_a
# mu = structured_singular_value(M)
# a = bisect_a(P, K, w; W = [], tol=1e-4)
# au = bisect_a(P, K, w; W = [], tol=1e-4, upper=true, N=2560)

# plot(w, mu, xscale=:log10, lab="mu")
# plot!(w, inv.(a), xscale=:log10, lab="lower")
# plot!(w, inv.(au), xscale=:log10, lab="upper")

## ==================================================


## 
# w = 2π .* exp10.(LinRange(-2, 2, 300))
# a = [-0.2 10;-10 -0.2]; b = I(2); c = [1 8;-10 1];
# P = ss(a,b,c,0);
# K = ss([1 -2;0 1]);
# dm = diskmargin(K*P) # disk margins at plant inputs
# dm = diskmargin(P*K); # disk margins at plant outputs
# MMIO = diskmargin(P,K,0,w)
# plot(MMIO.simultaneous, lab="simultaneous")
# plot!(MMIO.simultaneous_output, lab="output")
# plot!(MMIO.simultaneous_input, lab="input")

# ## Compare mu and diskmargin
# w = 2π .* exp10.(LinRange(-3, 3, 300))
# M,D = RobustAndOptimalControl.get_M(P, K, w; σ=0)
# mu = structured_singular_value(M)
# a = bisect_a(P, K, w; tol=1e-3, σ=0)
# au = bisect_a(P, K, w; tol=1e-3, upper=true, N=1024, σ=0)
# MMIO = diskmargin(P,K,0,w)

# dma = [dm.α for dm in MMIO.simultaneous]
# ##
# plot(w, (dma), lab="margin")
# plot!(w, inv.(mu), xscale=:log10, lab="mu", sp=1, l=:dash)
# plot!(w, a, xscale=:log10, lab="lower", sp=1)
# plot!(w, au, xscale=:log10, lab="upper", sp=1)

## Loop scaling

P = let
    tempA = [1.0 0.0 9.84 0.0 0.0; 0.0 1.0 0.01 2.14634e-6 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 -1.73983959887; 0.0 0.0 0.0 0.0 0.56597684805]
    tempB = [0.0 4.901416e-5 0.00019883877999999998; 0.0 0.0 0.0; 0.0 0.0 4.0414389999999996e-5; 0.0 -0.02004649371 0.0; 0.0 -0.00490141631 0.0]
    tempC = [0.0 -0.83516488404 0.0 0.0 0.0; 186.74725411661 0.0 0.0 0.0 0.0; -7.44299057498 0.0 7035.08410814126 0.0 0.0]
    ss(tempA, tempB, tempC, 0, 0.01)
end |> d2c
K, _ = glover_mcfarlane(P)

#

Li = K*P
Lo = P*K

So = output_sensitivity(P, K)

D = Diagonal([1000, 1, 1])
So = inv(D)*So*D
Sos = loop_scale(So, 30)

@test hinfnorm2(Sos)[1] <= 0.1*hinfnorm2(So)[1]

if isinteractive()
    w = exp10.(LinRange(-2, 2, 200))
    f1 = sigmaplot(So, w, c=1)
    sigmaplot!(Sos, w, c=2)

    f2 = muplot(So, w, c=1)
    muplot!(Sos, w, c=2)
    plot(f1,f2)
    display(current())
end




## Passivity

G = tf(1,[1,1]) |> ss
@test ispassive(G)
@test passivity_index(G) <= 1
passivityplot(G)

G = tf([1,1],[1,2]) |> ss
@test ispassive(G)
@test passivity_index(G) <= 1
passivityplot(G)

G = tf([1,2],[1,1]) |> ss
@test ispassive(G)
@test passivity_index(G) <= 1
passivityplot(G)


G = DemoSystems.resonant()
@test !ispassive(G)
@test passivity_index(G) > 1
passivityplot(G)


P1 = let
    P1A = [0.0 1.0; -1.0 -2.0]
    P1B = [0.0; 2.0;;]
    P1C = [-2.0 -3.5]
    P1D = [5.0;;]
    ss(P1A, P1B, P1C, P1D)
end
P2 = let
    P2A = [0.0 1.0 0.0; 0.0 0.0 2.0; -2.0 -1.5 -2.0]
    P2B = [0.0; 0.0; 2.0;;]
    P2C = [-0.975 0.5 -0.5]
    P2D = [1.0;;]
    ss(P2A, P2B, P2C, P2D)
end
P3 = P2*P1
@test ispassive(P1)
@test ispassive(P2)
@test !ispassive(P3)
passivityplot([P1, P2, P3])