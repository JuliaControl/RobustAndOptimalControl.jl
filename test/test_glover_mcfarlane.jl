using RobustAndOptimalControl, ControlSystems

G = tf(200, [10, 1])*tf(1, [0.05, 1])^2     |> ss
Gd = tf(100, [10, 1])                       |> ss
W1 = tf([1, 2], [1, 1e-6])                  |> ss
Gs = G*W1
Ks, γmin = glover_mcfarlane(Gs, 1.1)
@test γmin ≈ 2.34 atol=0.005

if isinteractive()
bodeplot([G, Gs, Gs*Ks]) |> display

plot( step(Gd*feedback(1, G*W1), 3))
plot!(step(Gd*feedback(1, G*W1*Ks), 3)) |> display

nyquistplot([G*W1, G*W1*Ks], ylims=(-2,1), xlims=(-2, 1), Ms_circles=1.5) |> display
end