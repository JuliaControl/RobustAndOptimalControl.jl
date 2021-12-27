using ControlSystems, RobustAndOptimalControl, MonteCarloMeasurements
using RobustAndOptimalControl: bisect_a

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
K = ss(1.0I(2))

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

dm = diskmargin(L3, ss(1.0I(3), L3.Ts), 0, w)
plot(dm.simultaneous)
plot!(dm.simultaneous_input)
plot!(dm.simultaneous_output)

plot(dm.input)
plot!(dm.output)

## MIMO

a = 10
P = ss([0 a; -a 0], I(2), [1 a; -a 1], 0)
K = ss(1.0I(2))

w = 2π .* exp10.(LinRange(-2, 2, 300))
# @time bisect_a(P, K, w)
##

# break at input (pass outputs through)
a = bisect_a(P, K, w; Z = [], tol=1e-4)
@test minimum(a) < 0.0998 # from DM paper

# break at output (pass inputs through)
a = bisect_a(P, K, w; W = [], tol=1e-4)
@test minimum(a) < 0.0998

# break at both input and output
a = bisect_a(P, K, w; tol=1e-4)
au = bisect_a(P, K, w; tol=1e-4, N=640, upper=true)
@test minimum(a) < 0.0499

ar = bisect_a(P, K, w; tol=1e-4, δ=δr)
aur = bisect_a(P, K, w; tol=1e-4, N=640, upper=true, δ=δr)

dm = diskmargin(P, K, 0, w)
adm = [dm.α for dm in dm.simultaneous]


if isinteractive()
    plot(w,  a,   xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Lower Complex")
    plot!(w, au,  xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Upper Complex")
    plot!(w, ar,  xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Lower Real")
    plot!(w, aur, xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="Upper Real")
    plot!(w, adm, xscale=:log10, xlabel="Frequency", ylims=(0,3), lab="μ")
    display(current())
end


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
M,_ = RobustAndOptimalControl.get_M(P, K, w; Z = :)


# Test stability of the computation
# @time mu = [structured_singular_value(M) for i = 1:10]
# plot(w, mu, xscale=:log10)

a = bisect_a(P, K, w; Z = [], tol=1e-4)
au = bisect_a(P, K, w; Z = [], tol=1e-4, upper=true, N=2560)
plot!(w, inv.(a), xscale=:log10)
plot!(w, inv.(au), xscale=:log10)

## ==================================================


## 
# NOTE: SISO och loop at a time blir rätt, men inte simultaneous. structured_singular_value verkar rätt, så är nog fel på feedback i get_M
w = 2π .* exp10.(LinRange(-2, 2, 300))
a = [-0.2 10;-10 -0.2]; b = I(2); c = [1 8;-10 1];
P = ss(a,b,c,0);
K = ss([1 -2;0 1]);
dm = diskmargin(K*P) # disk margins at plant inputs
dm = diskmargin(P*K); # disk margins at plant outputs
MMIO = diskmargin(P,K,0,w)

plot(MMIO.simultaneous)

w = 2π .* exp10.(LinRange(-3, 3, 300))

##
plot(MMIO.simultaneous, lab="simultaneous")
plot!(MMIO.simultaneous_output, lab="output") # denna är rätt för små frekvenser men 2x fel för höga
plot!(MMIO.simultaneous_input, lab="input") # samma för denna
