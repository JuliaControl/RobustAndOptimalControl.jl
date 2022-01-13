using RobustAndOptimalControl, ControlSystems

G = tf(200, [10, 1])*tf(1, [0.05, 1])^2     |> ss
Gd = tf(100, [10, 1])                       |> ss
W1 = tf([1, 2], [1, 1e-6])                  |> ss
Gs = G*W1
Ks, γ, info = glover_mcfarlane(Gs, 1.1)
@test info.γmin ≈ 2.34 atol=0.005

if isinteractive()
    bodeplot([G, Gs, Gs*Ks]) |> display

    plot( step(Gd*feedback(1, G*W1), 3))
    plot!(step(Gd*feedback(1, G*W1*Ks), 3)) |> display

    nyquistplot([G*W1, G*W1*Ks], ylims=(-2,1), xlims=(-2, 1), Ms_circles=1.5) |> display
end


W1h = hanus(W1)
@test W1h.nu == 2W1.nu

##

P = tf([1, 5], [1, 2, 10])
W1 = tf(1,[1, 0])

K,γ,info = glover_mcfarlane(ss(P); W1)


Gcl0 = let
    tempA = [-2.0 -2.5 2.0 0.0 0.0 0.0; 4.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -0.7189 0.0453 -2.9175; 0.0493 0.0616 0.0 -2.0493 -2.5616 2.0; 0.399 0.4988 0.0 3.601 -0.4988 0.0; 0.5 0.625 0.0 -1.2189 -0.5797 -2.9175]
    tempB = [0.0 2.0; 0.0 0.0; 0.0 0.0; 0.0985 0.0; 0.7981 0.0; 1.0 0.0]
    tempC = [0.5 0.625 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0]
    tempD = [1.0 0.0; 0.0 0.0]
    ss(tempA, tempB, tempC, tempD)
end

K0 = let
    tempA = [0.0 -0.7189 0.0453 -2.9175; 0.0 -2.0493 -2.5616 2.0; 0.0 3.601 -0.4988 0.0; 0.0 -1.2189 -0.5797 -2.9175]
    tempB = [0.0; 0.0985; 0.7981; 1.0;;]
    tempC = [1.0 0.0 0.0 0.0]
    tempD = [0.0;;]
    ss(tempA, tempB, tempC, tempD)
end

Gcl = extended_gangoffour(P, K)
@test hinfnorm2(Gcl-Gcl0)[1] < 1e-4
@test hinfnorm2(info.Gcl-Gcl0)[1] < 1e-4
@test nugap(K, -K0)[1] < 1e-4
@test info.margin ≈ 0.6325 atol=1e-3

@test ncfmargin(P, K)[1] ≈ 0.4472 atol=1e-3

p = ss(tf(4, [1, -0.001]))
cL = 1		
cH = 10
@test ncfmargin(p,cL)[1] ≈ 0.7069 atol=1e-3



P = ss(tf(1, [1,1]))
K = ss(1)
S, PS, CS, T = RobustAndOptimalControl.gangoffour2(P,K)
# gangoffourplot(P, K)
# bodeplot!(extended_gangoffour(P, K), plotphase=false)

isstable.((S, PS, CS, T))
gof = extended_gangoffour(P, K)
isstable(gof)
@test nugap(S, gof[1,1])[1] < 1e-6
@test nugap(PS, gof[1,2])[1] < 1e-6
@test nugap(CS, -gof[2,1])[1] < 1e-6 # NOTE: slightly disturbing to have - here
@test nugap(T, -gof[2,2])[1] < 1e-6
