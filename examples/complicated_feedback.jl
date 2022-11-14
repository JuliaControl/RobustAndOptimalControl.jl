#=
This example illustrates how named systems can be used to form complicated feedback interconnections.
=#
using RobustAndOptimalControl, ControlSystemsBase
const ROC = RobustAndOptimalControl
w = exp10.(LinRange(-2, 2, 300))
F = named_ss(ssrand(1, 1, 2, proper=true), x=:xF, u=:uF, y=:yF)
R = named_ss(ssrand(1, 1, 2, proper=true), x=:xR, u=:uR, y=:yR)
C = named_ss(ssrand(1, 1, 2, proper=true), x=:xC, u=:uC, y=:yC)
P = named_ss(ssrand(1, 1, 3, proper=true), x=:xP, u=:uP, y=:yP)

addP = sumblock("uP = yF + yC")
addC = sumblock("uC = yR - yP")

```
                 yF
              ┌────────────────────────────────┐
              │                                │
    ┌───────┐ │  ┌───────┐ yR   ┌─────────┐    │    ┌───────┐
uF  │       │ │  │       ├──────►         │ yC │  uP│       │    yP
────►   F   ├─┴──►   R   │      │    C    ├────+────►   P   ├────┬────►
    │       │    │       │   ┌──►         │         │       │    │
    └───────┘    └───────┘   │- └─────────┘         └───────┘    │
                             │                                   │
                             └───────────────────────────────────┘
```



# y1 = [:yP]
# u1 = [:uF]
# G = ROC.connect([F, R, C, P, addP, addC], u1, y1, z1=:)

connections = [
    :yP => :yP
    :uP => :uP
    :yC => :yC
    :yF => :yF
    :yF => :uR
    :uC => :uC
    :yR => :yR
]
w1 = [:uF]

G = ROC.connect([F, R, C, P, addP, addC], connections; w1)


@test sminreal(G[:yF, :uF].sys) ≈ F.sys
@test tf(sminreal(G[:yR, :uF].sys)) ≈ tf((R*F).sys)


# uF -> uP
manual = feedback(P.sys, C.sys)*F.sys + feedback(P.sys*C.sys)*R.sys*F.sys
automatic = G[:yP, :uF]
@test linfnorm(manual-automatic)[1] < 1e-6
isinteractive() && bodeplot([automatic, manual], w)
@test automatic.nx == 9

# uF -> yC
manual = feedback(C.sys, P.sys)*R.sys*F.sys - feedback(P.sys * C.sys)*F.sys
# manual = C*inv(1+P.sys*C.sys)*R.sys*F.sys
automatic = sminreal(G[:yC, :uF])
@test linfnorm(manual-automatic)[1] < 1e-6
isinteractive() && bodeplot([automatic, manual], w)
@test automatic.nx == 9



##

C = named_ss(ss([pid(2,1,form=:parallel) 0; 0 pid(5,6,form=:parallel)]), u=:e, y=:u, x=:xC)
G = named_ss(ss(-1,[1 2],[1;-1],0), x=:xG)

Sum = named_ss(ss([I(2) -I(2)]), u=[:r1, :r2, :y1, :y2], y=:e)
systems = [G,C,Sum]
u1 = [:r1, :r2]
y1 = [:y1, :y2]
T = ROC.connect(systems; u1, y1, w1=u1)


## manual
y1 = [:yP, :uP, :yC, :yF, :yF, :uC, :yR]
u1 = [:yP, :uP, :yC, :yF, :uR, :uC, :yR]
w1 = [:uF]

z2 = []
w2 = []

