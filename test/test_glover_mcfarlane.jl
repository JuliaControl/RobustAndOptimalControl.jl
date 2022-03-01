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


## Reduction
e,_ = ncfmargin(Gs,Ks)
Kr, hs, infor = baltrunc_coprime(Ks, n=Ks.nx)
n = findlast(RobustAndOptimalControl.error_bound(hs) .> 2e/3)
Ksr, hs, infor = baltrunc_coprime(Ks; n)
@test ncfmargin(Gs, Ksr)[1] >= 2/3 * e

## Discrete case
disc(G) = c2d(G, 0.01)
G = tf(200, [10, 1])*tf(1, [0.05, 1])^2     |> ss |> disc
Gd = tf(100, [10, 1])                       |> ss |> disc
W1 = tf([1, 2], [1, 1e-6])                  |> ss |> disc
K, γ, info = glover_mcfarlane(G, 1.1; W1)
@test info.γmin ≈ 2.4086 atol=0.005

Gcl0 = let
    tempA = [0.6542 -0.2067 -1.6744 0.2619 -0.3391 -1.0397 -6.283 -0.6716; 0.1309 0.9823 -0.1434 0.0224 -0.029 -0.089 -0.538 -0.0575; 0.0014 0.0199 0.999 0.0002 -0.0002 -0.0006 -0.0037 -0.0004; 0.0 0.0 -0.1271 1.0 -0.0259 -0.0794 -0.4799 -0.0513; 0.0 0.0 -1.8079 0.0 0.315 -1.2464 -6.1494 -0.4097; 0.0 0.0 0.3374 0.0 0.1019 0.8933 -1.0187 -0.0351; 0.0 0.0 0.1454 0.0 0.0012 0.0193 0.8499 -0.0002; 0.0 0.0 0.0181 0.0 -0.0259 -0.0794 -0.6251 0.9487]
    tempB = [-0.1065 0.1309; -0.0091 0.0112; -0.0001 0.0001; -0.0081 0.0; -0.1157 0.0; 0.0216 0.0; 0.0093 0.0; 0.0012 0.0]
    tempC = [0.0 0.0 15.625 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -12.7105 2.0 -2.5903 -7.941 -47.9883 -5.1294]
    tempD = [1.0 0.0; -0.8135 0.0]
    ss(tempA, tempB, tempC, tempD, 0.01)
end

@test hinfnorm2(Gcl0 - info.Gcl)[1] < 0.02

if isinteractive()
    bodeplot([G, G*K]) |> display
    plot( step(Gd*feedback(1, info.Gs), 3))
    plot!(step(Gd*feedback(1, G*K), 3)) |> display
    nyquistplot([info.Gs, G*K], ylims=(-2,1), xlims=(-2, 1), Ms_circles=1.5) |> display
end

Ko = observer_controller(info)
@test nugap(G*K, info.Gs*Ko)[1] < 1e-6 # The loop transfer is used since the observer controller is for the scaled plant whereas K is for the unscaled plant.
# bodeplot([G*K, info.Gs*Ko])


## Test the observer predictor
Gs = info.Gs
Ph = observer_predictor(info)
u = randn(Gs.nu, 100)
res = lsim(Gs, u)
uy = [u; res.y]
reso = lsim(Ph, uy)


# plot(res)
# plot!(reso) |> display
# plot(res.x')
# plot!(reso.x') |> display

@test res.y ≈ reso.y
@test res.x ≈ reso.x

## Strictly proper case
Ksp, γsp, infosp = glover_mcfarlane(G, 1.1; W1, strictly_proper=true)
@test infosp.γmin > 2.4086
@test iszero(Ksp.D)

Ko = observer_controller(infosp)
@test nugap(G*Ksp, infosp.Gs*Ko)[1] < 1e-6
# bodeplot([G*Ksp, infosp.Gs*Ko])

# Can be manually verified to be "slightly worse" than the controller without strictly proper restriction.
# if isinteractive()
#     plot( step(Gd*feedback(1, info.Gs), 3))
#     plot!(step(Gd*feedback(1, G*K), 3))
#     plot!(step(Gd*feedback(1, G*Ksp), 3)) |> display
#     nyquistplot([info.Gs, G*K, G*Ksp], ylims=(-2,1), xlims=(-2, 1), Ms_circles=1.5) |> display
# end

Gs = info.Gs
Ph = observer_predictor(infosp)
u = randn(Gs.nu, 100)
res = lsim(Gs, u)
uy = [u; res.y]
reso = lsim(Ph, uy)
@test res.y ≈ reso.y
@test res.x ≈ reso.x


##
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

@test all(isstable.((S, PS, CS, T)))
gof = extended_gangoffour(P, K)
isstable(gof)
@test nugap(S, gof[1,1])[1] < 1e-6
@test nugap(PS, gof[1,2])[1] < 1e-6
@test nugap(CS, -gof[2,1])[1] < 1e-6 # NOTE: slightly disturbing to have - here
@test nugap(T, -gof[2,2])[1] < 1e-6


## 2 DOF GMcF

P = tf([1, 5], [1, 2, 10])
W1 = tf(1,[1, 0]) |> ss

Tref = tf(1, [1, 1]) |> ss

K1dof, γ1, info1 = glover_mcfarlane(ss(P), 1.1; W1)
K2dof, γ2, info2 = RobustAndOptimalControl.glover_mcfarlane_2dof(ss(P), Tref, 1.1, 1.1; W1)

G1 = feedback(P*K1dof)
G2 = info2.Gcl

plot([step(G1, 15), step(G2, 15), step(Tref, 15)], lab=["1-DOF" "2-DOF" "Tref"])
w = exp10.(LinRange(-3, 2, 200))
bodeplot(info2.K1, w, lab="Feedforward filter")

@test dcgain(G2)[] ≈ 1 rtol=1e-4