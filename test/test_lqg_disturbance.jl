using RobustAndOptimalControl, ControlSystemsBase, LinearAlgebra, Test

@testset "Observer-only constant disturbances" begin
    for discrete in (false, true)
        P = ss(tf(1, [1, 1, 0]))
        discrete && (P = c2d(P, 0.2))
        Pa = add_low_frequency_disturbance(P)
        pe = ExtendedStateSpace(Pa; B1=I, C1=[P.C zeros(P.ny, P.nu)])
        prob = LQGProblem(pe, [2.0], [0.3], ones(Pa.nx), [0.1];
            disturbance_indices=[Pa.nx], disturbance_input_indices=[1], qQ=0.2, qR=0.1)
        L, K = lqr(prob), kalman(prob)
        Lx = lqr(P, 2.2P.C'P.C, [0.3;;])
        @test L[:, 1:P.nx] ≈ Lx
        @test L[:, end:end] == ones(1, 1)
        @test system_mapping(prob, identity).nx == P.nx
        @test observer_predictor(prob).nx == Pa.nx
        @test observer_filter(prob).nx == Pa.nx
        Cfb = observer_controller(prob)
        manual = observer_controller(Pa, [Lx ones(1, 1)], K)
        @test freqresp(Cfb, [0.1, 1, 10]) ≈ freqresp(manual, [0.1, 1, 10])
        @test isstable(feedback(P, Cfb))
        @test norm(dcgain(G_PS(prob))) < 1e-8
        @test norm(step(G_PS(prob), 100).y[:, end, :]) < 1e-6
        @test isstable(closedloop(prob))
        @test dcgain(closedloop(prob)*static_gain_compensation(prob)) ≈ ones(1, 1)
        Cff = ff_controller(prob)
        @test dcgain(minreal(feedback(P, Cfb)*Cff, 1e-7)) ≈ ones(1, 1)
        @test freqresp(closedloop(prob), [0.1, 1, 10]) ≈ freqresp(feedback(P, Cfb)*ff_controller(prob; comp_dc=false), [0.1, 1, 10])
        for (actual, expected) in zip(gangoffour(prob), gangoffour(P, Cfb))
            @test freqresp(actual, [0.1, 1, 10]) ≈ freqresp(expected, [0.1, 1, 10])
        end
        Ce, cl = extended_controller(prob; z=[1])
        @test isstable(cl)
        @test freqresp(system_mapping(Ce), [0.1, 1]) ≈ -freqresp(Cfb, [0.1, 1])
        @test size(lft(prob)) == (1, Pa.nx)
        @test freqresp(lft(prob), [0.1, 1]) ≈ freqresp(lft(pe, -Cfb), [0.1, 1])
        @test_throws ArgumentError lqi_controller(prob, ones(1, 1))
        @test performance_mapping(prob).nu == Pa.nx
        inferred = LQGProblem(pe, [2.0], [0.3], ones(Pa.nx), [0.1]; disturbance_indices=[Pa.nx], qQ=0.2)
        @test lqr(inferred) ≈ L
        if discrete
            Cd = observer_controller(prob; direct=true)
            @test isstable(feedback(P, Cd))
            @test norm(dcgain(G_PS(P, Cd))) < 1e-8
        end
    end

    @testset "MIMO ordering and explicit allocation" begin
        P = ss([-1.0 0; 0 -2], [1.0 0; 0 2], Matrix{Float64}(I, 2, 2), 0)
        E = [0.0 1; 1 0]
        Pa = add_low_frequency_disturbance(P, P.B*E)
        perm = [3, 1, 4, 2]
        Pa = ss(Pa.A[perm, perm], Pa.B[perm, :], Pa.C[:, perm], Pa.D)
        pe = ExtendedStateSpace(Pa; B1=I)
        SQ = zeros(4, 2); SQ[[2, 4], :] = 0.1I(2)
        prob = LQGProblem(pe, [1.0, 2], [0.2, 0.4], ones(4), ones(2);
            disturbance_indices=[1, 3], disturbance_input_indices=[2, 1], SQ, qQ=0.3)
        @test lqr(prob)[:, [1, 3]] == E
        @test lqr(prob)[:, [2, 4]] ≈ lqr(P, Diagonal([1.3, 2.3]), Diagonal([0.2, 0.4]), SQ[[2, 4], :])
        @test isstable(feedback(P, observer_controller(prob)))
        @test norm(dcgain(G_PS(prob))) < 1e-8
        @test dcgain(closedloop(prob)*static_gain_compensation(prob)) ≈ I(2)

        for kwargs in ((disturbance_indices=[1, 1],), (disturbance_indices=[5],),
                       (disturbance_indices=[1, 3], disturbance_input_indices=[1]),
                       (disturbance_indices=[1, 3], disturbance_input_indices=[1, 3]),
                       (disturbance_indices=[1, 3], disturbance_input_indices=[1, 2]))
            @test_throws ArgumentError LQGProblem(pe, ones(2), ones(2), ones(4), ones(2); kwargs...)
        end
        A = copy(pe.A); A[1, 1] = -0.001
        @test_throws ArgumentError LQGProblem(ExtendedStateSpace(pe; A), ones(2), ones(2), ones(4), ones(2); disturbance_indices=[1, 3])
        C1 = copy(pe.C1); C1[1, 1] = 1
        @test_throws ArgumentError LQGProblem(ExtendedStateSpace(pe; C1), ones(2), ones(2), ones(4), ones(2); disturbance_indices=[1, 3])


        for overrides in (
            (; A=pe.A + [1.0; 0; 0; 0] * [0.0 1 0 0]),
            (; B2=pe.B2 + [1.0; 0; 0; 0] * [1.0 0]),
            (; C2=pe.C2 + [1.0; 0] * [1.0 0 0 0]),
        )
            bad = ExtendedStateSpace(pe; overrides...)
            @test_throws ArgumentError LQGProblem(bad, ones(2), ones(2), ones(4), ones(2); disturbance_indices=[1, 3])
        end
        @test_throws ArgumentError LQGProblem(pe, ones(2), ones(2), ones(4), ones(2); disturbance_indices=[1, 3], SQ=ones(4, 2))
        @test_throws ArgumentError LQGProblem(pe, ones(2), ones(2), ones(4), ones(2); disturbance_input_indices=[1])
        unmatched = ss([-1.0 0; 0 -2], [1.0; 0;;], Matrix{Float64}(I, 2, 2), 0)
        unmatched = add_low_frequency_disturbance(unmatched, [0.0; 1;;])
        @test_throws ArgumentError LQGProblem(ExtendedStateSpace(unmatched; B1=I), ones(2), ones(1), ones(3), ones(2); disturbance_indices=[3])

        redundant = ss([-1.0;;], [1.0 1], [1.0;;], 0)
        Pa = add_low_frequency_disturbance(redundant, redundant.B[:, 2:2])
        pe = ExtendedStateSpace(Pa; B1=I)
        @test_throws ArgumentError LQGProblem(pe, ones(1), ones(2), ones(2), ones(1); disturbance_indices=[2])
        prob = LQGProblem(pe, ones(1), ones(2), ones(2), ones(1); disturbance_indices=[2], disturbance_input_indices=[2])
        @test lqr(prob)[:, 2] == [0, 1]
    end
end
