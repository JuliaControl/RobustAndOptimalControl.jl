using RobustAndOptimalControl, ControlSystems

# example 19.1
A = [
    -1 0 4
    0 2 0
    0 0 -3
]
B1 = [
    0 0
    1 0
    0 0
]
C1 = B1'
B2 = [
    1
    1
    1
]
C2 = B2'
D11 = 0I(2)
D12 = [
    0
    1
]
D21 = D12'
D22 = 0I(1)

P = ExtendedStateSpace(A, B1, B2, C1, C2, D11, D12, D21, D22)
for γ ∈ (nothing, 1000)
    K, Cl = h2synthesize(P, γ) # the example uses h2syn
    @test norm(Cl) ≈ 8.9582 rtol = 1e-3

    K_ML = let
        A = [-1.0 -4.236148004700622 4.0; -4.23607198996497 -6.472215758593602 -4.23607198996497; 0.0 -4.236148004700622 -3.0]
        B = [0.0; 4.23607198996497; 0.0]
        C = [0.0 -4.236148004700622 0.0]
        D = [0.0]
        sys = ss(A,B,C,D)
    end
    @test K ≈ K_ML rtol=1e-3



    Cl_ML = let
        A = [-1.0 0.0 4.0 0.0 -4.236148004700622 0.0; 0.0 2.0 0.0 0.0 -4.236148004700622 0.0; 0.0 0.0 -3.0 0.0 -4.236148004700622 0.0; 0.0 0.0 0.0 -1.0 -4.236148004700622 4.0; 4.23607198996497 4.23607198996497 4.23607198996497 -4.23607198996497 -6.472215758593602 -4.23607198996497; 0.0 0.0 0.0 0.0 -4.236148004700622 -3.0]
        B = [0.0 0.0; 1.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 4.23607198996497; 0.0 0.0]
        C = [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -4.236148004700622 0.0]
        D = [0.0 0.0; 0.0 0.0]
        ss(A,B,C,D)
    end

    @test Cl_ML ≈ Cl rtol=1e-3

    # s = tf("s")
    # K⁺ = ss(-148.79(s+1)*(s+3) / ((s+31.74)*(s+3.85)))
    # RobustAndOptimalControl.controller_reduction_weight(P, K⁺)
    # K̂⁺ = controller_reduction(P, K⁺, 1, true)
    # @test K̂⁺ ≈ tf(-117.085, [1, 34.526])
end


# Test conversion between H2 problem and LQG
K, Cl = h2synthesize(P)
#=
    The controller above is given by -L(sI-Ae)\K where L is feedback gain and K is Kalman gain.
    We have the following equivalences

    C1'C1    = Q
    D12'D12  = R

    B1*B1'   = R1
    D21*D21' = R2
=#
Q = C1'C1
R = D12'D12
R1 = B1*B1'
R2 = D21*D21'

L = lqr(system_mapping(P), Q, R)
@test L ≈ -K.C

Kal = kalman(system_mapping(P), R1, R2)
@test Kal ≈ K.B

l = LQGProblem(P)
@test lqr(l) == L
@test kalman(l) == Kal

# Implement normalized coprime fact from sec 13.8