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

    K_matlab = let
        A = [-1.0 -4.236148004700622 4.0; -4.23607198996497 -6.472215758593602 -4.23607198996497; 0.0 -4.236148004700622 -3.0]
        B = [0.0; 4.23607198996497; 0.0]
        C = [0.0 -4.236148004700622 0.0]
        D = [0.0]
        sys = ss(A,B,C,D)
    end
    @test K ≈ K_matlab rtol=1e-3



    Cl_matlab = let
        A = [-1.0 0.0 4.0 0.0 -4.236148004700622 0.0; 0.0 2.0 0.0 0.0 -4.236148004700622 0.0; 0.0 0.0 -3.0 0.0 -4.236148004700622 0.0; 0.0 0.0 0.0 -1.0 -4.236148004700622 4.0; 4.23607198996497 4.23607198996497 4.23607198996497 -4.23607198996497 -6.472215758593602 -4.23607198996497; 0.0 0.0 0.0 0.0 -4.236148004700622 -3.0]
        B = [0.0 0.0; 1.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 4.23607198996497; 0.0 0.0]
        C = [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -4.236148004700622 0.0]
        D = [0.0 0.0; 0.0 0.0]
        ss(A,B,C,D)
    end

    @test Cl_matlab ≈ Cl rtol=1e-3

    # s = tf("s")
    # K⁺ = ss(-148.79(s+1)*(s+3) / ((s+31.74)*(s+3.85)))
    # RobustAndOptimalControl.controller_reduction_weight(P, K⁺)
    # K̂⁺ = controller_reduction(P, K⁺, 1, true)
    # @test K̂⁺ ≈ tf(-117.085, [1, 34.526])
end


# TODO: implement LQG in terms of ExtendedStateSpace. C1 = M, B1 = N

# Q = C1'C1
# R = I with u = sqrt(R₀)u₀
# Use this go from weighted specification to standard Q/R which works with KalmanFilter? Verify by comparing Riccati equation 14.9 to that of the lqr function

# Implement normalized coprime fact from sec 13.8