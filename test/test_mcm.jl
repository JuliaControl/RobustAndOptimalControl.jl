using RobustAndOptimalControl
import RobustAndOptimalControl as RC
using MonteCarloMeasurements
using ControlSystems
unsafe_comparisons(false)
G = tf(1, [1 ± 0.1, 1])



funs = [
    G->c2d(G, 0.01)
    pole
    tzero
    minreal
    ss
    bode
    step
    impulse
    dcgain
]

for fun in funs
    @show fun
    @test_nowarn fun(G)
end 

@test dcgain(G)[] == 1
@test real(pole(G)[]) == -denvec(G)[][1]


H = ss(G)

funs = [
    G->c2d(G, 0.01)
    pole
    tzero
    minreal
    balreal
    tf
    # bode
    # step
    # impulse
    dcgain
]


for fun in funs
    @show fun
    @test_nowarn fun(H)
end 

RC.sys2Rn(pole, H)


w = exp10.(LinRange(-1, 1, 100))
@test_nowarn bode(H, w)