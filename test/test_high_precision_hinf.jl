using GenericLinearAlgebra

## Test higher prec hinfsyn
# Example from "A robust numerical method for the γ-iteration in H∞ control"
a = b = d = n = e2 = 1.0
e1 = 0
A = -1.0I(2)
B1 = [e1 0; 0 e2]
B2 = ones(2)
C1 = diagm([a, b])
C2 = [d n]
D11 = diagm([0.5, 0.5])
D12 = [0; 1]
D21 = [0 1]
D22 = [0;;]
P = ss(A,B1,B2,C1,C2,D11,D12,D21,D22)
K, γ = hinfsynthesize(P, γrel=1.01, ftype=BigFloat)[1:2] # 
@test_broken hinfnorm2(lft(P, K)) ≈ 0.806 atol=1e-2
@test_broken γ ≈ 0.806 atol=1e-2 # Not numerically robust enough to pass this test despite highprec
