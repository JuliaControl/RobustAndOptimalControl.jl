using RobustAndOptimalControl
using ControlSystemsBase
using Test
using LinearAlgebra
using Random


using RobustAndOptimalControl: rqr

@testset "rqr" begin
    @info "Testing rqr"

    e = sqrt(eps())
    D = [1 1; e 0; 0 e]
    γ = 0.00001
    b = randn(2)

    bD = big.(D)
    bb = big.(b)
    hp = (bD'bD + big(γ)*I)\bb
    
    @test norm(rqr(D, γ)\b - hp) < 1e-10
    @test norm(rqr(D, γ)\b - hp) < norm((D'D + γ*I)\b - hp)

    @test rqr(D, γ)*b ≈ (D'D + γ*I)*b


    γ = 0.00001
    b = randn(3,2)
    bb = big.(b)
    hp = bb/(bD'bD + big(γ)*I)
    
    @test norm(b/rqr(D, γ) - hp) < 1e-10
    @test norm(b/rqr(D, γ) - hp) < norm(b/(D'D + γ*I) - hp)

    @test b*rqr(D, γ) ≈ b*(D'D + γ*I)

    display(rqr(D, γ))
end


"""
Tests for the public and private methods of the hInfSynthesis function. This
function utilizes the preexisting ControlSystemsBase toolbox, and performs a
H-infinity synthesis using the dual Riccati equation approach. As such,
the synthesis is done in a set of steps.

(1) Re-writing the specifications on an extended state-space form.
(2) Verifying that the resulting extended state-space object satisfies a set of
assumptions required for proceeding with the synthesis.
(3) A coordinate transform to enable the synthesis.
(4) Synthesis using the γ-iterations, checking if a solution to the H-infinity
problem exists in each iteration and applying a bisection method.
(5) Re-transforming the system to the original coordinates
(6) Verification that the computed solution is correct.

In addition to these six ponts, the code also enables

(7) A bilinear discretization with an inverse operation to move from continuous
to discrete time, thereby enabling approximate discrete-time synthesis.
(8) Plotting functionality to visualize the H-infinity synthesis.
(9) Three examples which can be used to demonstrate the tool.

Many of the tests are quite intuitive, and all points (1)-(9) are tested
extensively with detailed comments for each test-set, and below is a list
of all the functions tested in this unit test

hinfsynthesize
hinfassumptions
hinfsignals
hinfpartition
bilineard2c
bilinearc2d
_detectable
_stabilizable
_synthesizecontroller
_assertrealandpsd
_checkfeasibility
_solvehamiltonianare
_solvematrixequations
_γIterations
_transformp2pbar
_scalematrix
_input2ss
_coordinateTtansformqr
_coordinateTtansformsvd

"""

@testset "(1) Specifications" begin
    """
    These tests make sure that the specifications are written on a suitable form
    for any combination of none/empty, static gain, SISO, and MIMO
    specificaitons in the weighting functions. Essentially, the tests are made
    to verify that hInf_partition(G, WS, WU, WT) works as expected.
    """

    @testset "Conversion of user input" begin

        @testset "Empty and nothing" begin

            # Check that conversion of nothing is OK
            A, B, C, D = RobustAndOptimalControl._input2ss(nothing)
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test size(A) == (0, 0)
            @test size(B) == (0, 0)
            @test size(C) == (0, 0)
            @test size(D) == (0, 0)

            # Check that conversion of empty objects are OK
            A, B, C, D = RobustAndOptimalControl._input2ss([])
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test size(A) == (0, 0)
            @test size(B) == (0, 0)
            @test size(C) == (0, 0)
            @test size(D) == (0, 0)
        end

        @testset "Static gains" begin
            # Fixture
            number = 2.0

            # Check that conversion of numbers are OK
            A, B, C, D = RobustAndOptimalControl._input2ss(number)
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test D[1, 1] == number

            @test size(A) == (0, 0)
            @test size(B) == (0, 1)
            @test size(C) == (1, 0)
            @test size(D) == (1, 1)
        end

        @testset "LTI models" begin

            # Fixture
            Random.seed!(0)
            M = 3
            N = 2
            SISO_tf = tf([1, 0], [1, 0, 1])
            SISO_ss = ss(SISO_tf)
            MIMO_ss = ss(
                rand(Float64, (M, M)),
                rand(Float64, (M, N)),
                rand(Float64, (N, M)),
                rand(Float64, (N, N)),
            )
            MIMO_tf = tf(MIMO_ss)

            # check that conversion of SISO tf data is OK
            A, B, C, D = RobustAndOptimalControl._input2ss(SISO_tf)
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test ss(A, B, C, D) == ss(SISO_tf)

            # check that conversion of SISO tf data is OK
            A, B, C, D = RobustAndOptimalControl._input2ss(SISO_ss)
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test ss(A, B, C, D) == SISO_ss

            # check that conversion of MIMO tf data is OK
            A, B, C, D = RobustAndOptimalControl._input2ss(MIMO_tf)
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test ss(A, B, C, D) == ss(MIMO_tf)

            # check that conversion of MIMO tf data is OK
            A, B, C, D = RobustAndOptimalControl._input2ss(MIMO_ss)
            @test isa(A, Array{Float64,2})
            @test isa(B, Array{Float64,2})
            @test isa(C, Array{Float64,2})
            @test isa(D, Array{Float64,2})

            @test ss(A, B, C, D) == MIMO_ss
        end
    end

    #
    # @testset "Dimensionality checks"
    #
    #   # Check sensitivity function weight - must have the same number of inputs
    #   # as outputs, can have an arbitrarily large statespace, and must have the
    #   # same number of inputs as there are outputs in the process which is to be
    #   # controlled.
    #
    #   # Fixture
    #   M = 5; N = 3; L=4;
    #   G = ss(rand(Float64, (M,M)), rand(Float64, (M,L)),
    #          rand(Float64, (N,M)), rand(Float64, (N,L)))
    #
    #   # Create a sensitivity weighting function with ii inputs and ii outputs
    #   for ii = 2:6
    #     for jj = 2:6
    #       for kk = 2:6
    #         # Create some randome weighting functions
    #         WS = ss(rand(Float64, (M,M)), rand(Float64, (M,ii)),
    #                 rand(Float64, (ii,M)), rand(Float64, (ii,ii)))
    #         WU = ss(rand(Float64, (M,M)), rand(Float64, (M,jj)),
    #                 rand(Float64, (jj,M)), rand(Float64, (jj,jj)))
    #         WT = ss(rand(Float64, (M,M)), rand(Float64, (M,kk)),
    #                 rand(Float64, (kk,M)), rand(Float64, (kk,kk)))
    #
    #         println([ii,jj,kk])
    #         # Chech that the specifications can be re-written is possible
    #         if ii == N && jj == L && kk == N
    #           @test isa(hInf_partition(G, WS, WU, WT), ControlSystemsBase.ExtendedStateSpace)
    #         else
    #           @test_throws ErrorException hInf_partition(G, WS, WU, WT)
    #         end
    #       end
    #     end
    #   end
    # end
end


@testset "(2) Assumptions" begin
    """
    Tests the methods used to check that the assumptions are satisfied,
    incorporating two separate algorithms to determine detectability and
    observability, and a rectangular pseudoinverse to check the rank conditions.
    """

    @testset "Stabilizability check" begin
        """
        Test the check for stabilizability using the Hautus Lemma
        """
        ### Fixture
        Random.seed!(0)
        N = 10
        M = 5
        R = rand(Float64, (N, N))
        Q = eigvecs(R + R')
        L = rand(Float64, N) .- 0.5
        A = Q * Diagonal(L) * Q'
        E = eigvals(A)

        # Here we should return false if for any positive eigenvalue λ of A, we get
        # that rank(A-λI, B) < N. We define B as conlums which are linearly
        # dependent with columns in A-λI, and check that the system returns false
        # in every case where λ > 0 and never when λ≦0.
        for ii = 1:N
            for jj = 1:(N-M)
                Ahat = A - E[ii] * Matrix{Float64}(I, N, N)
                B = Ahat[:, jj:(jj+M)]
                if E[ii] >= 0
                    @test !RobustAndOptimalControl._stabilizable(A, B)
                else
                    @test RobustAndOptimalControl._stabilizable(A, B)
                end
            end
        end

        # Check common input types and ensure that a method error is thrown when
        # not using two abstract matrices, but some other common type
        N = 10
        M = 5
        A = rand(Float64, (N, N))
        B = rand(Float64, (N, M))
        @test_throws MethodError RobustAndOptimalControl._stabilizable(A, nothing)
        @test_throws MethodError RobustAndOptimalControl._stabilizable(nothing, B)
        @test_throws MethodError RobustAndOptimalControl._stabilizable(A, [])
        @test_throws MethodError RobustAndOptimalControl._stabilizable([], B)
        @test_throws MethodError RobustAndOptimalControl._stabilizable(A, ss(1))
        @test_throws MethodError RobustAndOptimalControl._stabilizable(ss(1), B)
        @test_throws MethodError RobustAndOptimalControl._stabilizable(A, tf(1))
        @test_throws MethodError RobustAndOptimalControl._stabilizable(tf(1), B)
    end

    @testset "Detectability check" begin
        """
        Test the check for detectability using the Hautus Lemma
        """
        ### Fixture
        Random.seed!(0)
        N = 10
        M = 5
        R = rand(Float64, (N, N))
        Q = eigvecs(R + R')
        L = rand(Float64, N) .- 0.5
        A = Q * Diagonal(L) * Q'
        E = eigvals(A)

        # Here we should return false if for any positive eigenvalue λ of A, we get
        # that rank(A-λI; C) < N. We define C as rows which are linearly
        # dependent with rows in A-λI, and check that the system returns false
        # in every case where λ > 0 and never when λ≦0.
        for ii = 1:N
            for jj = 1:(N-M)
                Ahat = A - E[ii] * Matrix{Float64}(I, N, N)
                C = Ahat[jj:(jj+M), :]
                if E[ii] >= 0
                    @test !RobustAndOptimalControl._detectable(A, C)
                else
                    @test RobustAndOptimalControl._detectable(A, C)
                end
            end
        end

        # Check common input types and ensure that a method error is thrown when
        # not using two abstract matrices, but some other common type
        N = 10
        M = 5
        A = rand(Float64, (M, M))
        C = rand(Float64, (N, M))
        @test_throws MethodError RobustAndOptimalControl._detectable(A, nothing)
        @test_throws MethodError RobustAndOptimalControl._detectable(nothing, C)
        @test_throws MethodError RobustAndOptimalControl._detectable(A, [])
        @test_throws MethodError RobustAndOptimalControl._detectable([], C)
        @test_throws MethodError RobustAndOptimalControl._detectable(A, ss(1))
        @test_throws MethodError RobustAndOptimalControl._detectable(ss(1), C)
        @test_throws MethodError RobustAndOptimalControl._detectable(A, tf(1))
        @test_throws MethodError RobustAndOptimalControl._detectable(tf(1), C)
    end

    # TODO: write tests using the above submethods directly in hInf_assumptions
end


@testset "(3) Coordinate transform" begin
    """
    Computes the coordinate transfomration to write the system on the form of
    assumption A4 in [1], i.e. take an extended statespace model, P, satisfying
    assumption A2, and transform is to a system P_bar, where the matrices
    D12_bar = [0,I] and D21_bar = [0;I].
    """

    @testset "Find transformation for full rank matrix with QR" begin
        """
        Test computation of the transformation in square and rectangular cases.
        """
        # Fixture
        tolerance = 1e-10

        for N = 1:5
            for M = 1:5
                A = rand(Float64, (M, N))
                I_mat = Matrix{Float64}(I, min(M, N), min(M, N))

                if M == N
                    # Square case
                    A_bar_true = I_mat
                elseif M > N
                    # Rectangular case (more rows than columns)
                    Z_mat = zeros(Float64, (max(M, N) - min(M, N), min(M, N)))
                    A_bar_true = [Z_mat; I_mat]
                else
                    # Rectangular case (more columns than rows)
                    Z_mat = zeros(Float64, (min(M, N), max(M, N) - min(M, N)))
                    A_bar_true = [Z_mat I_mat]
                end
                Tl, Tr = RobustAndOptimalControl._coordinatetransformqr(A)
                @test opnorm(Tl * A * Tr - A_bar_true) < tolerance
            end
        end
    end

    @testset "Find transformation for full rank matrix with SVD" begin
        """
        Test computation of the transformation in square and rectangular cases.
        """
        # Fixture
        tolerance = 1e-10

        for N = 1:5
            for M = 1:5
                A = rand(Float64, (M, N))
                I_mat = Matrix{Float64}(I, min(M, N), min(M, N))

                if M == N
                    # Square case
                    A_bar_true = I_mat
                elseif M > N
                    # Rectangular case (more rows than columns)
                    Z_mat = zeros(Float64, (max(M, N) - min(M, N), min(M, N)))
                    A_bar_true = [Z_mat; I_mat]
                else
                    # Rectangular case (more columns than rows)
                    Z_mat = zeros(Float64, (min(M, N), max(M, N) - min(M, N)))
                    A_bar_true = [Z_mat I_mat]
                end
                Tl, Tr = RobustAndOptimalControl._coordinatetransformsvd(A)
                @test opnorm(Tl * A * Tr - A_bar_true) < tolerance
            end
        end
    end

    @testset "Find transformation for arbitrary matrix inputs" begin
        """
        Test computation of the transformation in square and rectangular cases.
        """
        # Fixture
        test_matrices_full = [
            rand(Float64, (1, 1)),
            rand(Float64, (5, 5)),
            rand(Float64, (2, 5)),
            rand(Float64, (5, 2)),
        ]

        # Compute rank deficient matrices
        test_matrices_rank_deficient = test_matrices_full
        for ii = 1:length(test_matrices_full)
            U, S, V = svd(test_matrices_full[ii])
            S[1] = 0
            test_matrices_rank_deficient[ii] = U * Diagonal(S) * V'
        end

        # Make sure that an exception is raised if assumption A2 is violated,
        # any of the matrices which are to be decomposed are rank deficient
        for A in test_matrices_rank_deficient
            println(A)
            @test_throws ErrorException RobustAndOptimalControl._scalematrix(A)
        end

        # Various bad inputs
        @test_throws MethodError RobustAndOptimalControl._scalematrix(tf(1))
        @test_throws MethodError RobustAndOptimalControl._scalematrix(ss(1))
        @test_throws MethodError RobustAndOptimalControl._scalematrix([])
        @test_throws MethodError RobustAndOptimalControl._scalematrix(nothing)
        @test_throws MethodError RobustAndOptimalControl._scalematrix(1)

        # Check that errors are thrown if the matric has zise zero
        @test_throws ErrorException RobustAndOptimalControl._scalematrix(zeros(Float64, (0, 1)))

        # Various bad methods
        A = test_matrices_full[1]
        @test_throws ErrorException RobustAndOptimalControl._scalematrix(
            A;
            method = "bad method",
        )
        @test_throws ErrorException RobustAndOptimalControl._scalematrix(A; method = 1)
    end

    @testset "Test application of the transformation to an ESS object" begin
        """
        The wa
        """
    end
end


@testset "(4) Gamma iterations" begin
    """
    Tests the core methods of the γ-iteration bisection method, including
    the ARE hamiltonian solver, the eigenvalue solver, the feasibility checks.
    """

    @testset "Solution feasibility check" begin
        """
        Check that a solution to the dual riccati equations is correctly reported
        as being feasible if satisfying the conditions of positive definiteness on
        the X and Y solutions are met, and the spectral radius condition of X*Y is
        also met
        """
        # Fixture
        tolerance = 1e-10
        iteration = 1
        Random.seed!(0)
        N = 10
        M = 5
        R = rand(Float64, N, N)
        Q = eigvecs(R + R')
        ρX = 0.2
        ρY = 0.3
        LX = rand(Float64, N)
        LX = ρX * sort(LX / maximum(LX))
        LY = rand(Float64, N)
        LY = ρY * sort(LY / maximum(LY))
        Xinf = Q * Diagonal(LX) * Q'
        Yinf = Q * Diagonal(LY) * Q'
        γ = 1

        # Test that the fesibility is true if ρ(Xinf*Yinf) < γ^2, that is the
        # check should be true for any
        #
        #    γ = sqrt(ρX*ρY) + ϵ
        #
        # with any epsilon greater than or equal to zero.
        @test !RobustAndOptimalControl._checkfeasibility(
            Xinf,
            Yinf,
            sqrt(ρX * ρY) - tolerance,
            tolerance,
            iteration;
            verbose = false,
        )
        @test RobustAndOptimalControl._checkfeasibility(
            Xinf,
            Yinf,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )

        # Test that errors are thrown if the matrix Xinf and Yinf are not PSD down
        # to the numerical tolerance.
        L = LX
        L[1] += -L[1] + 2 * tolerance
        Xpos = Q * Diagonal(L) * Q' # slightly positive eigenvalue
        L[1] += -L[1]
        Xzero = Q * Diagonal(LX) * Q'               # exactly one zero eigenvalue
        L[1] += -L[1] - 2 * tolerance
        Xneg = Q * Diagonal(L) * Q'   # slightly negative eigenvalue
        @test RobustAndOptimalControl._checkfeasibility(
            Xpos,
            Yinf,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )
        @test RobustAndOptimalControl._checkfeasibility(
            Xzero,
            Yinf,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )
        @test !RobustAndOptimalControl._checkfeasibility(
            Xneg,
            Yinf,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )

        L = LY
        L[1] += -L[1] + 2 * tolerance
        Ypos = Q * Diagonal(L) * Q' # slightly positive eigenvalue
        L[1] += -L[1]
        Yzero = Q * Diagonal(L) * Q'               # exactly one zero eigenvalue
        L[1] += -L[1] - 2 * tolerance
        Yneg = Q * Diagonal(L) * Q'   # slightly negative eigenvalue
        @test RobustAndOptimalControl._checkfeasibility(
            Xinf,
            Ypos,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )
        @test RobustAndOptimalControl._checkfeasibility(
            Xinf,
            Yzero,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )
        @test !RobustAndOptimalControl._checkfeasibility(
            Xinf,
            Yneg,
            sqrt(ρX * ρY) + tolerance,
            tolerance,
            iteration;
            verbose = false,
        )
    end

    # TODO: Include a check to verify that the bisection works as intended.
    # TODO: Check to verify that the Hamiltonian Shur-solver is working
    # TODO: Check to verify that the Hamiltonian Eigenvalue-solver is working
end


@testset "(7) Bilinear discretization" begin
    """
    This tests the bilinear method of discretizing a continuous time StateSpace
    or ExtendedStateSpace object, moving from the Laplace-domain to the Z-domain.
    However, importantly, the method also enables approximating discrete-time
    systems with continuous-time state-space domain, moving from the Z-domain back
    to the laplace-domain. Since the L∞-norm is invariant over the bilinear
    discretization and "continuization", this effectively allows us to approximate
    discrete-time systems with a continuous-time equivalent, design an H∞-optimal
    controller for the continuous-time plant, and then re-discretizing this
    controller.

    The structure of these tests is to simply discretize a set of plants with
    various dimensions, transport them back into continuous-time domain, and then
    to discrete-time again. By then comparing the resulting two pairs of discrete-
    time and continuous-time systems in the frequency domain, we can test if the
    bilinear discretization operates as expected.
    """

    @testset "SS data" begin
        # Fixture for the SS-type data
        Random.seed!(0)
        N = [1, 5, 10]
        M = [1, 5, 10]
        H = [0.1, 0.01, 0.001]
        tolerance = 1e-7

        # Tests for the SS-type data
        for (ii, m) in enumerate(N)
            for (jj, n) in enumerate(M)
                for (kk, h) in enumerate(H)
                    freq = [10^i for i in range(-6, stop = log10(pi / h), length = 101)]

                    AcTrue, BcTrue, CcTrue, DcTrue =
                        rand(m, m), rand(m, n), rand(n, m), rand(n, n)

                    ev = eigvals(AcTrue)
                    if !isempty(ev[real(ev).<0])
                        AcTrue =
                            AcTrue -
                            one(AcTrue) * 2 * maximum(abs.(real(ev[real(ev).<0])))
                    end

                    Gtrue = ss(AcTrue, BcTrue, CcTrue, DcTrue)
                    valA, fA = sigma(Gtrue, freq)
                    sysB = bilinearc2d(Gtrue, h)
                    valB, fB = sigma(sysB, freq)
                    sysC = bilineard2c(bilinearc2d(Gtrue, h))
                    valC, fC = sigma(sysC, freq)
                    sysD = bilinearc2d(bilineard2c(bilinearc2d(Gtrue, h)), h)
                    valD, fD = sigma(sysD, freq)

                    @test abs(maximum(svd(sysB.B).S) - maximum(svd(sysB.C).S)) <
                            tolerance
                    @test abs(maximum(svd(sysC.B).S) - maximum(svd(sysC.C).S)) <
                            tolerance
                    @test abs(maximum(svd(sysD.B).S) - maximum(svd(sysD.C).S)) <
                            tolerance

                    # Test that the C->D->C and D->C->D results in the same
                    @test norm(valA - valC, Inf) < tolerance
                    @test norm(valB - valD, Inf) < tolerance
                end
            end
        end
    end

    @testset "Extended SS data" begin
        # Fixture for the extended SS-type data
        Random.seed!(0)
        N = [1, 4]
        P1 = [1, 4]
        P2 = [1, 4]
        M1 = [1, 4]
        M2 = [1, 4]
        H = [0.1, 0.01, 0.001]
        tolerance = 1e-5

        for (in1, n) in enumerate(N)
            for (im1, m1) in enumerate(M1)
                for (im2, m2) in enumerate(M2)
                    for (ip1, p1) in enumerate(P1)
                        for (ip2, p2) in enumerate(P2)
                            for (ih, h) in enumerate(H)
                                freq = [
                                    10^i for i in
                                    range(-4, stop = log10(pi / h), length = 101)
                                ]

                                A, B1, B2 = 0.01 * rand(n, n), rand(n, m1), rand(n, m2)
                                C1, D11, D12 = rand(p1, n), rand(p1, m1), rand(p1, m2)
                                C2, D21, D22 = rand(p2, n), rand(p2, m1), rand(p2, m2)

                                ev = eigvals(A)
                                if !isempty(ev[real(ev).>0])
                                    A =
                                        A -
                                        one(A) *
                                        2 *
                                        maximum(abs.(real(ev[real(ev).>0])))
                                end

                                EsysA = ss(A, B1, B2, C1, C2, D11, D12, D21, D22)
                                
                                valA, fA = sigma(
                                    ss(
                                        EsysA.A,
                                        [EsysA.B1 EsysA.B2],
                                        [EsysA.C1; EsysA.C2],
                                        [EsysA.D11 EsysA.D12; EsysA.D21 EsysA.D22],
                                    ),
                                    freq,
                                )

                                # To discrete time
                                EsysB = bilinearc2d(EsysA, h)
                                valB, fB = sigma(
                                    ss(
                                        EsysB.A,
                                        [EsysB.B1 EsysB.B2],
                                        [EsysB.C1; EsysB.C2],
                                        [EsysB.D11 EsysB.D12; EsysB.D21 EsysB.D22],
                                        h,
                                    ),
                                    freq,
                                )

                                # To continuous time
                                EsysC = bilineard2c(bilinearc2d(EsysA, h))
                                valC, fC = sigma(
                                    ss(
                                        EsysC.A,
                                        [EsysC.B1 EsysC.B2],
                                        [EsysC.C1; EsysC.C2],
                                        [EsysC.D11 EsysC.D12; EsysC.D21 EsysC.D22],
                                    ),
                                    freq,
                                )

                                # To discrete time
                                EsysD =
                                    bilinearc2d(bilineard2c(bilinearc2d(EsysA, h)), h)
                                valD, fD = sigma(
                                    ss(
                                        EsysD.A,
                                        [EsysD.B1 EsysD.B2],
                                        [EsysD.C1; EsysD.C2],
                                        [EsysD.D11 EsysD.D12; EsysD.D21 EsysD.D22],
                                        h,
                                    ),
                                    freq,
                                )

                                # Test that the C->D->C and D->C->D results in the same
                                @test norm(valA - valC, Inf) < tolerance
                                @test norm(valB - valD, Inf) < tolerance
                            end
                        end
                    end
                end
            end
        end
    end
end



@testset "(9) Synthesis Examples" begin

    @testset "DC motor example" begin
        # Fixture
        tolerance = 1e-2

        # Make sure that the code runs
        include("../examples/hinf_example_DC.jl")
        Ω = [10^i for i in range(-3, stop = 3, length = 201)]

        # Check that the optimal gain is correct
        @test_broken abs(γ - 2.1982) < tolerance

        # Check that the closed loop satisfies ||F_l(P(jω), C(jω)||_∞ < γ,  ∀ω ∈ Ω
        valPcl = sigma(Pcl, Ω)[1]
        @test all(valPcl .< (γ + tolerance))

        # Check that ||S(jω)/WS(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WS, LTISystem) || isa(WS, Number)
            valSWS = sigma(S * WS, Ω)[1]
            @test all(valSWS .< (γ + tolerance))
        end

        # Check that ||C(jω)S(jω)/W_U(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WU, LTISystem) || isa(WU, Number)
            valKSWU = sigma(CS * WU, Ω)[1]
            @test all(valKSWU .< (γ + tolerance))
        end

        # Check that ||T(jω)/W_T(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WT, LTISystem) || isa(WT, Number)
            valTWT = sigma(T * WT, Ω)[1]
            @test all(valTWT .< (γ + tolerance))
        end

        # C_test = let
        #     _A = [-1.99970880438311e-7 -3.28465977959578e-12 -1.46200072004469e-10; 5213.86459190004 -599.125731164477 -26175.7449027424; 0.000132853900209503 0.249985014202232 -0.000667020485484034]
        #     _B = [3.9999999999998; -5.48401518999795e-14; -9.03638143773611e-7]
        #     _C = [651.733073987504 -74.8757163955595 -3271.9681128428]
        #     _D = [0.0]
        #     ss(_A, _B, _C, _D)
        # end
        # @test C ≈ C_test

        # Pcl_test = let
        #     _A = [-2.0e-7 0.0 -22.4 0.0 0.0 0.0; 0.0 -0.12 0.0 5213.86459190003 -599.005731164476 -26175.7449027424; 0.0 0.25 0.0 0.0 0.0 0.0; 0.0 0.0 -22.3999999999989 -1.99970880438311e-7 -3.28465977959578e-12 -1.46200072004469e-10; 0.0 0.0 3.07104850639885e-13 5213.86459190004 -599.125731164477 -26175.7449027424; 0.0 0.0 5.06037360513222e-6 0.000132853900209503 0.249985014202232 -0.000667020485484034]
        #     _B = [4.0; 0.0; 0.0; 3.9999999999998; -5.48401518999795e-14; -9.03638143773611e-7]
        #     _C = [4.99999996666667 0.0 -3.73333333333333 0.0 0.0 0.0; 0.0 0.0 0.0 651.733073987504 -74.8757163955595 -3271.9681128428]
        #     _D = [0.666666666666667; 0.0]
        #     ss(_A, _B, _C, _D)
        # end
        # @test Pcl ≈ Pcl_test


        include("../examples/hinf_example_DC_discrete.jl")
    end

    @testset "MIT open courseware example" begin

        # Fixture
        tolerance = 1e-3

        # Make sure that the code runs
        @test_nowarn include("../examples/hinf_example_MIT.jl")
        Ω = [10^i for i in range(-7, stop = 7, length = 201)]

        # Check that the optimal gain is correct
        @test abs(γ - 0.923430124918) < tolerance

        # Check that the closed loop satisfies ||F_l(P(jω), C(jω)||_∞ < γ,  ∀ω ∈ Ω
        valPcl = sigma(Pcl, Ω)[1]
        @test all(valPcl .< (γ + tolerance))

        # Check that ||S(jω)/WS(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WS, LTISystem) || isa(WS, Number)
            valSWS = sigma(S * WS, Ω)[1]
            @test all(valSWS .< (γ + tolerance))
        end

        # Check that ||C(jω)S(jω)/W_U(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WU, LTISystem) || isa(WU, Number)
            valKSWU = sigma(CS * WU, Ω)[1]
            @test all(valKSWU .< (γ + tolerance))
        end

        # Check that ||T(jω)/W_T(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WT, LTISystem) || isa(WT, Number)
            valTWT = sigma(T * WT, Ω)[1]
            @test all(valTWT .< (γ + tolerance))
        end
        @test_nowarn include("../examples/hinf_example_MIT_discrete.jl")
    end

    @testset "Quad tank example" begin

        # Fixture
        tolerance = 1e-2

        # Make sure that the code runs
        @test_nowarn include("../examples/hinf_example_tank.jl")
        Ω = [10^i for i in range(-7, stop = 7, length = 201)]

        # Check that the optimal gain is correct
        @test abs(γ - 0.9534467) < tolerance

        # Check that the closed loop satisfies ||F_l(P(jω), C(jω)||_∞ < γ,  ∀ω ∈ Ω
        valPcl = sigma(Pcl, Ω)[1]
        @test all(valPcl[:, 1] .< (γ + tolerance))

        # Check that ||S(jω)/WS(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WS, LTISystem) || isa(WS, Number)
            valSWS = sigma(S * WS, Ω)[1]
            @test all(valSWS[:, 1] .< (γ + tolerance))
        end

        # Check that ||C(jω)S(jω)/W_U(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WU, LTISystem) || isa(WU, Number)
            valKSWU = sigma(CS * WU, Ω)[1]
            @test all(valKSWU[:, 1] .< (γ + tolerance))
        end

        # Check that ||T(jω)/W_T(jω)||_∞ < γ,  ∀ω ∈ Ω
        if isa(WT, LTISystem) || isa(WT, Number)
            valTWT = sigma(T * WT, Ω)[1]
            @test all(valTWT[:, 1] .< (γ + tolerance))
        end
        @test_nowarn include("../examples/hinf_example_tank_discrete.jl")
    end
end



# Active suspension model
Ghinf = let
    tempA = [-0.0002624650388438681 -0.0007086326853279555 -5.514860214838735e-6 5.392822831615305e-6 -0.00021941859259636337 5.7367174759374515e-5 2.716271789984803e-6 -5.207362655859199e-8; 0.0007085461243516112 -0.0015661701446945551 -3.6743299740329606e-5 4.446857259890745e-5 -0.0012415813360395955 0.00031572689863442214 1.4919863256506027e-5 3.774700671281847e-7; -2.6688550308198786e-5 0.00011969041716615464 -0.005470740837813468 -3.312837378399592 -0.0006648000498561453 7.736877960108244e-5 8.048249950431425e-5 -0.0007036720083433111; 2.6780262946570124e-5 -0.00012072390764621185 3.3128423564838365 -0.005488882784039457 -5.4875578252832274e-5 -8.643382344753049e-5 7.21668048337983e-5 0.0006989710246561228; -0.00021786003811079174 0.0012362483346265432 0.000518056606377057 -0.0002817704182696309 -0.010171435814673755 0.005284664019198872 0.00025699262090504574 0.00017835246885565745; -5.714887667033384e-6 2.4538056618037784e-5 -0.0023358039105166465 0.0023437022418646533 -8.220448017148232e-5 -0.014400492240553908 3.3120788303644457 -0.485883016448181; -3.43900972024107e-7 1.4955212699039632e-6 -0.0001372374389515997 0.0001350438735793037 -3.1648355641875196e-6 -3.313608890780306 -4.876599282157503e-5 -0.025115291089941788; -0.000606695465422999 0.0026543480323023747 -0.24781110296757258 0.24863240120600844 0.020904866034982067 -2.9780981841413716 -0.1738223493687497 -994.9864755102109]
    tempB = [-2.6916440270002493 -0.2554009951371452; 2.92286785112518 0.22008351310595206; 0.0001087824228703122 -1.4605771184212994; 0.0003553184094320467 1.4605799810748008; -1.1361512541642613 0.04015266842029391; 7.351988871983848e-5 -0.31162757907120153; -5.5192817054244374e-5 -0.018134024020639343; -1.4581188719236578e-5 -33.08054048313464]
    tempC = [-1.6662263861541294e-5 7.489883766321597e-5 -1.0325586314158979 -1.0325495438963512 -0.00018999604391594334 -0.0003438875147835997 0.0072442285059155925 -5.8212476684194005e-8; -2.703733193978162 -2.931128701591983 -0.0295545993317383 0.030476380583964593 -1.13522384486694 0.29547249708088114 0.013984829519988138 -0.00014480521475524215; 2.0171842012718553e-5 -8.82532858129867e-5 0.008266731048337915 -0.008239080945429876 -0.0006950217342690853 0.09903294293255775 0.005449445704509105 33.08054048282091; -0.0020128336425470905 0.008817549966645846 -1.032553335277771 -1.0325399529992152 -0.06097714335767885 -0.00033892067271467107 0.007148084470341972 -6.909583102265332e-8]
    tempD = [0.0 0.0; 0.0 0.0; 0.0 1.1; 0.001 0.0]
    named_ss(ss(tempA, tempB, tempC, tempD), x=[:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8], u=[:do, :u], y=[Symbol("err₊y(t)_scaled"), :e_scaled, :uw_scaled, :y])
end

Gsyn = partition(Ghinf, u = [:u], y = [:y])
hinfassumptions(Gsyn)
K, γ, mats = hinfsynthesize(Gsyn, γrel = 1.02, transform=true, verbose=true, method = :SVD)
@test isstable(K)
@test γ < 10.000042479470356 * 1.01

# Gcl, S, CS, T = hinfsignals(Gsyn, Ghinf[:err, :u], K)

# using Plots
# gangoffourplot(Ghinf[:err, :u].sys, K)