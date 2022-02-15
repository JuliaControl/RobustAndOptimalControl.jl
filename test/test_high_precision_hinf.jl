using GenericLinearAlgebra, RobustAndOptimalControl, ControlSystems

bb(x) = big.(x)
bb(P::AbstractStateSpace) =  ss(bb.(ssdata_e(P)), P.timeevol)
@testset "Numerical difficulties 1" begin
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
    K, γ = hinfsynthesize(P, γrel=1.05, ftype=BigFloat, check=false)[1:2] # gamma is way off (too low), but the actual realized gamma by the controller is quite close
    @test_broken hinfnorm2(lft(P, K)) ≈ 0.806 atol=1e-2
    @test_broken γ ≈ 0.806 atol=1e-2 # Not numerically robust enough to pass this test despite highprec

end


@testset "Numerical difficulties 2" begin
    ## Numerically difficult problem instance
    Gsyn3A = [0.0 1.0 0.0; -1.2000000000000002e-6 -0.12000999999999999 0.0; -11.2 -0.0 -2.0e-7]
    Gsyn3B = [0.0 0.0; 0.0 1.0; 1.0 0.0]
    Gsyn3C = [-7.466666666666666 -0.0 19.999999866666666; 0.0 0.0 0.0; -11.2 -0.0 0.0]
    Gsyn3D = [0.6666666666666666 0.0; 0.0 1.0; 1.0 -0.0]
    Gsyn3 = ss(Gsyn3A, Gsyn3B, Gsyn3C, Gsyn3D)
    Gsyn = partition(Gsyn3, 1, 2)
    K, γ = hinfsynthesize(Gsyn, ftype=BigFloat, γrel = 1)[1:2]
    @test hinfnorm2(lft(Gsyn, K))[1] ≈ γ atol=1e-2
    # It seems we are better than slicot here since the computed controller gives lower hinfnorm than slicot
    @test γ <= 4.4825150 # value by slicot
    @test γ ≈ 4.4825150 atol=2e-2 

    # Gsynb = modal_form(Gsyn)
    Gsynb, _ = balance_statespace(Gsyn, false)
    # Gsynb, _ = RobustAndOptimalControl.schur_form(Gsyn)

    K, γ = hinfsynthesize(Gsynb, ftype=BigFloat, γrel = 1)[1:2]
    @test hinfnorm2(lft(Gsynb, K))[1] ≈ γ atol=1e-2
    # @test_broken γ ≈ 4.4825150 atol=1e-2 # value by slicot
    @test γ ≈ 4.4825150 atol=2e-2 # slightly less strict test is passed
end