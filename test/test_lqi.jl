using Test
using ControlSystemsBase
using RobustAndOptimalControl
using LinearAlgebra

"""
    designed_poles(G, K, L, inds; ϵ = 0)

The closed-loop poles that an LQI design is supposed to realize: by the separation principle they
are the poles of the augmented state-feedback loop together with the observer poles. Comparing
these against the poles of the system assembled by `lqi_controller` verifies not only the sign but
also the scaling of every channel of `L`, including the integrator channels.
"""
function designed_poles(G, K, L, inds; ϵ = 0)
    Ga = add_output_integrator(G, inds; ϵ)
    sortpoles([eigvals(Ga.A - Ga.B * L); eigvals(G.A - K * G.C)])
end
sortpoles(p) = sort(p, by = x -> (real(x), imag(x)))

"Close the loop around an `lqi_controller`, which already contains the negative feedback sign."
function lqi_loop(C, G; name = "plant")
    Gn = G isa NamedStateSpace ? G : named_ss(G, name = name, x = :x_plant, y = :y_plant, u = :u_plant)
    refs = [r for r in C.u if endswith(string(r), "_r")]
    feedback(C, Gn, w1 = refs, z2 = Gn.y, u1 = Gn.y, pos_feedback = true)
end


@testset "SISO continuous" begin
    G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
    Q = diagm([0,5])
    R = [1.0;;]
    K = kalman(G,Q,R)
    obs = observer_predictor(G,K; output_state=true)

    Q1 = diagm([0.488,0,100])
    Q2 = [1/100;;]
    L = lqi(G,Q1,Q2)

    # Pin the sign and the scaling of the gain against an explicitly augmented plant. The
    # integrator state integrates +y, so the augmented A has +C in the lower-left block.
    A, B, C, D = ssdata(G)
    Aa = [A zeros(2,1); C 0]
    Ba = [B; D]
    @test L ≈ lqr(ss(Aa, Ba, [C 0], D), Q1, Q2)

    C0 = lqi_controller(G, obs, Q1, Q2)

    @test C0.nu == 2
    @test C0.ny == G.nu
    @test count(p -> abs(p) < 1e-8, poles(C0)) == 1  # one integral mode

    H = lqi_loop(C0, G)
    @test sortpoles(poles(H)) ≈ designed_poles(G, K, L, [1])
    @test dcgain(H)[2] ≈ 1

    res = step(H, 50)
    @test res.y[:, end] ≈ [dcgain(feedback(-C0[:, 2], G))[]; 1.0]

    # The integrator rejects a static load disturbance entering at the plant input exactly
    @test dcgain(feedback(G, -ss(C0[:, 2])))[] ≈ 0 atol=1e-10
    @test dcgain(output_sensitivity(G, -ss(C0[:, 2])))[] ≈ 0 atol=1e-10
end

@testset "SISO continuous, ϵ > 0" begin
    G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
    K = kalman(G, diagm([0,5]), [1.0;;])
    obs = observer_predictor(G, K; output_state=true)
    Q1 = diagm([0.488,0,100]); Q2 = [1/100;;]
    ϵ = 0.5
    L = lqi(G, Q1, Q2; ϵ)
    A, B, C, D = ssdata(G)
    @test L ≈ lqr(ss([A zeros(2,1); C -ϵ], [B; D], [C 0], D), Q1, Q2)

    Cc = lqi_controller(G, obs, Q1, Q2; ϵ)
    H = lqi_loop(Cc, G)
    @test sortpoles(poles(H)) ≈ designed_poles(G, K, L, [1]; ϵ)
    # A finite integrator pole trades exact tracking for a bounded low-frequency controller gain
    @test !(dcgain(H)[2] ≈ 1)
    @test any(p -> p ≈ -ϵ, poles(Cc))
end

@testset "MIMO continuous, all outputs integrated" begin
    G2 = ss([-1.0 0.2; 0.1 -2.0], [1.0 0; 0 1.0], [1.0 0; 0 1.0], 0)
    K2 = kalman(G2, I(G2.nx), I(G2.ny))
    obs2 = observer_predictor(G2, K2; output_state=true)
    Q1_mimo = diagm([1.0, 1.0, 10.0, 10.0])
    Q2_mimo = Matrix{Float64}(I(G2.nu))
    L = lqi(G2, Q1_mimo, Q2_mimo)
    C_mimo = lqi_controller(G2, obs2, Q1_mimo, Q2_mimo)
    @test C_mimo.nu == 2*G2.ny       # [r; y]
    @test C_mimo.ny == G2.nu
    # Each output should have an integral mode in the controller
    @test count(p -> abs(p) < 1e-6, poles(C_mimo)) == G2.ny

    H2 = lqi_loop(C_mimo, G2)
    @test isstable(minreal(H2))
    @test sortpoles(poles(H2)) ≈ designed_poles(G2, K2, L, 1:2)
    # Reference-to-output DC gain should be identity on the plant outputs.
    H2_dc = dcgain(H2)
    y_out_inds = [findfirst(==(s), H2.y) for s in [:y_plant1, :y_plant2]]
    @test H2_dc[y_out_inds, :] ≈ I(G2.ny) atol=1e-8
end

@testset "Partial integrator_outputs" begin
    # integrate only output 1 on a 2-output plant
    G2 = ss([-1.0 0.2; 0.1 -2.0], [1.0 0; 0 1.0], [1.0 0; 0 1.0], 0)
    K2 = kalman(G2, I(G2.nx), I(G2.ny))
    obs2 = observer_predictor(G2, K2; output_state=true)
    Q2_mimo = Matrix{Float64}(I(G2.nu))
    Q1_partial = diagm([1.0, 1.0, 10.0])  # nx + 1 integrator
    L = lqi(G2, Q1_partial, Q2_mimo; integrator_outputs=[1])
    C_partial = lqi_controller(G2, obs2, Q1_partial, Q2_mimo; integrator_outputs=[1])
    # Controller takes [r_for_integrated_output; all y] = [1 ref; ny measurements]
    @test C_partial.nu == 1 + G2.ny
    @test count(p -> abs(p) < 1e-6, poles(C_partial)) == 1
    H_partial = lqi_loop(C_partial, G2)
    @test isstable(minreal(H_partial))
    @test sortpoles(poles(H_partial)) ≈ designed_poles(G2, K2, L, [1])
    H_partial_dc = dcgain(H_partial)
    y_out_inds_p = [findfirst(==(s), H_partial.y) for s in [:y_plant1, :y_plant2]]
    # Only the integrated output should track exactly
    @test H_partial_dc[y_out_inds_p[1], 1] ≈ 1 atol=1e-8

    # A scalar index is accepted and equivalent to the length-one vector
    @test lqi(G2, Q1_partial, Q2_mimo; integrator_outputs=1) ≈ L
    @test ss(lqi_controller(G2, obs2, Q1_partial, Q2_mimo; integrator_outputs=1)) ≈ ss(C_partial)
end

@testset "integrator_outputs order" begin
    # The integrator states, the reference channels and the integrator entries of Q1 must all
    # follow the order in which the indices are given. An asymmetric plant with asymmetric
    # integrator weights makes a swapped pairing detectable.
    G = ss([-1.0 0.0; 0.0 -2.0], [1.0 0; 0 1.0], [1.0 0; 0 3.0], 0)
    K = kalman(G, Matrix(1.0I,2,2), Matrix(1.0I,2,2))
    obs = observer_predictor(G, K; output_state=true)
    Q1 = diagm([1.0, 1.0, 1.0, 100.0])
    Q2 = Matrix(1.0I, 2, 2)

    for inds in ([1,2], [2,1])
        L = lqi(G, Q1, Q2; integrator_outputs=inds)
        Ga = add_output_integrator(G, inds)
        # Integrator state k integrates output inds[k]
        @test Ga.A[3:4, 1:2] ≈ G.C[inds, :]

        C = lqi_controller(G, obs, Q1, Q2; integrator_outputs=inds)
        @test C.u[1:2] == Symbol.("y_plant" .* string.(inds) .* "_r")
        H = lqi_loop(C, G)
        @test isstable(H)
        @test sortpoles(poles(H)) ≈ designed_poles(G, K, L, inds)
        # Each reference drives its own output to unit gain, whatever order it was given in
        dc = dcgain(H)
        yi = [findfirst(==(s), H.y) for s in [:y_plant1, :y_plant2]]
        ri = [findfirst(==(Symbol("y_plant$(i)_r")), H.u) for i in 1:2]
        @test dc[yi, ri] ≈ I(2) atol=1e-8
    end

    # The two orders correspond to different problems, since Q1 is read in the given order
    @test !(lqi(G, Q1, Q2; integrator_outputs=[1,2]) ≈ lqi(G, Q1, Q2; integrator_outputs=[2,1]))
end

@testset "Nonzero D" begin
    G = ss([-1.0 0.5; 0 -2], [0.0; 1.0;;], [1.0 0.0], 0.7)
    K = kalman(G, Matrix(1.0I,2,2), [1.0;;])
    obs = observer_predictor(G, K; output_state=true)
    Q1 = diagm([1.0, 1, 10]); Q2 = [1.0;;]
    L = lqi(G, Q1, Q2)
    # A nonzero D feeds the control signal into the integrator, so B_aug = [B; D]
    A, B, C, D = ssdata(G)
    @test L ≈ lqr(ss([A zeros(2,1); C 0], [B; D], [C 0], D), Q1, Q2)
    H = lqi_loop(lqi_controller(G, obs, Q1, Q2), G)
    @test sortpoles(poles(H)) ≈ designed_poles(G, K, L, [1])
    @test dcgain(H)[2] ≈ 1
end

@testset "Discrete" begin
    G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
    Ts = 0.1
    Gd = c2d(G, Ts)
    Kd = kalman(Gd, Matrix{Float64}(diagm([0,5])), [1.0;;])
    obsd = observer_predictor(Gd, Kd; output_state=true)
    Q1d = diagm([0.488, 0, 100.0])
    Q2d = [1/100;;]
    A, B, C, D = ssdata(Gd)

    for ϵd in (0.0, 1e-4)
        Ld = lqi(Gd, Q1d, Q2d; ϵ=ϵd)
        # The discrete integrator state is a forward-Euler time integral, xᵢ⁺ = (1-ϵ)xᵢ + Ts*y,
        # so the Ts factor belongs in the augmented A and B rather than in the added output.
        Aa = [A zeros(2,1); Ts*C (1-ϵd)]
        Ba = [B; Ts*D]
        @test Ld ≈ lqr(ss(Aa, Ba, [C 0], D, Ts), Q1d, Q2d)

        Cd = lqi_controller(Gd, obsd, Q1d, Q2d; ϵ=ϵd)
        @test Cd.nu == 2
        # Discrete integrator pole near 1 (offset by ϵ)
        @test any(p -> abs(p - (1 - ϵd)) < 1e-8, poles(Cd))
        Hd = lqi_loop(Cd, Gd)
        @test isstable(minreal(Hd))
        # The realized loop must be the designed loop, which fails if the integrator state of the
        # controller differs from that of the plant augmentation by a factor Ts
        @test sortpoles(poles(Hd)) ≈ designed_poles(Gd, Kd, Ld, [1]; ϵ=ϵd)
    end

    Cd = lqi_controller(Gd, obsd, Q1d, Q2d)
    Hd = lqi_loop(Cd, Gd)
    @test dcgain(Hd)[2] ≈ 1 atol=1e-6
    @test dcgain(feedback(Gd, -ss(Cd[:, 2])))[] ≈ 0 atol=1e-10
end

@testset "Cross term" begin
    G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
    K = kalman(G, diagm([0,5]), [1.0;;])
    obs = observer_predictor(G, K; output_state=true)
    Q1 = diagm([0.488, 0, 100.0]); Q2 = [1/100;;]
    S = 0.01*ones(3, 1)
    L = lqi(G, Q1, Q2, S)
    @test L ≈ lqr(add_output_integrator(G, [1]), Q1, Q2, S)
    @test !(L ≈ lqi(G, Q1, Q2))
    H = lqi_loop(lqi_controller(G, obs, Q1, Q2, S), G)
    @test sortpoles(poles(H)) ≈ designed_poles(G, K, L, [1])
end

@testset "LQGProblem method" begin
    G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
    Q = diagm([0,5]); R = [1.0;;]
    Q1 = diagm([0.488, 0, 100]); Q2 = [1/100;;]
    K = kalman(G, Q, R)
    obs = observer_predictor(G, K; output_state=true)
    C0 = lqi_controller(G, obs, Q1, Q2)

    prob = LQGProblem(G, diagm([0.488, 0]), Q2, Matrix{Float64}(Q), Matrix{Float64}(R))
    Qi = [100.0;;]
    C_prob = lqi_controller(prob, Qi)
    @test C_prob.nu == 2
    @test count(p -> abs(p) < 1e-8, poles(C_prob)) == 1
    # With C1 = I and qQ = 0 the LQGProblem method reproduces the explicit-observer form exactly
    @test C_prob.nx == C0.nx
    @test C_prob.ny == C0.ny
    @test ss(C_prob) ≈ ss(C0)
    @test_throws ArgumentError lqi_controller(prob, [100.0 0; 0 1.0])

    # `prob.Q1` penalizes the performance output C1*x, and qQ/SQ must be honoured, exactly as in
    # `lqr(::LQGProblem)`
    Pe = ExtendedStateSpace(G, B1=I(2), C1=[1.0 0.0; 0.0 2.0])
    for qQ in (0.0, 10.0)
        probe = LQGProblem(Pe, diagm([1.0, 3.0]), Q2, Matrix(1.0I,2,2), R; qQ)
        Ge = system_mapping(probe, identity)
        Ke = kalman(probe)
        Q1_aug = cat(probe.C1'probe.Q1*probe.C1 + qQ*probe.C2'probe.C2, Qi; dims=(1,2))
        Le = lqi(Ge, Q1_aug, probe.Q2, [probe.SQ; zeros(1, Ge.nu)])
        Ce = lqi_controller(probe, Qi)
        @test sortpoles(poles(lqi_loop(Ce, Ge))) ≈ designed_poles(Ge, Ke, Le, [1])
    end
    # A nonzero qQ changes the design rather than being silently discarded
    prob_q0 = LQGProblem(Pe, diagm([1.0, 3.0]), Q2, Matrix(1.0I,2,2), R; qQ=0.0)
    prob_q1 = LQGProblem(Pe, diagm([1.0, 3.0]), Q2, Matrix(1.0I,2,2), R; qQ=10.0)
    @test !(ss(lqi_controller(prob_q0, Qi)) ≈ ss(lqi_controller(prob_q1, Qi)))
end

@testset "Validation" begin
    G = ss([0 32;-31.25 -0.4],[0; 2.236068],[0.0698771 0],0)
    G2 = ss([-1.0 0.2; 0.1 -2.0], [1.0 0; 0 1.0], [1.0 0; 0 1.0], 0)
    Q1 = diagm([0.488, 0, 100.0]); Q2 = [1/100;;]

    @test_throws ArgumentError lqi(G, ones(3,4), Q2)                       # Q1 not square
    @test_throws ArgumentError lqi(G, diagm([1.0,1,1,1]), Q2)              # Q1 wrong size
    @test_throws ArgumentError lqi(G, Q1, ones(1,2))                       # Q2 not square
    @test_throws ArgumentError lqi(G, Q1, Q2; integrator_outputs=Int[])    # no integrator
    @test_throws ArgumentError lqi(G, Q1, Q2; integrator_outputs=[2])      # out of range
    @test_throws ArgumentError lqi(G2, diagm([1.0,1,1,1]), Matrix(1.0I,2,2); integrator_outputs=[1,1]) # duplicate

    # Integrating more outputs than there are control inputs leaves the augmentation
    # unstabilizable, which must be rejected rather than silently returning a useless gain
    Gwide = ss([-1.0 0; 0 -2], reshape([1.0, 1.0], 2, 1), [1.0 0; 0 1.0], 0)
    @test_throws ArgumentError lqi(Gwide, diagm([1.0,1,10,10]), [1.0;;])

    # A transmission zero at the integrator pole also destroys stabilizability, warn about it
    Gz = ss(tf([1.0, 0.0], [1.0, 3, 2]))
    @test_logs (:warn,) match_mode=:any try
        lqi(Gz, diagm([1.0, 1, 10]), [1.0;;])
    catch
    end

    # The observer must map [u; y] to the full state estimate
    K2 = kalman(G2, I(G2.nx), I(G2.ny))
    @test_throws ArgumentError lqi_controller(G2, observer_predictor(G2, K2; output_state=true)[1, :],
                                              diagm([1.0,1,10,10]), Matrix(1.0I,2,2))
end
