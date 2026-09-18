using RobustAndOptimalControl, ControlSystemsBase
G1 = ssrand(2,3,4, proper=true)
G2 = ssrand(2,3,4, proper=true)


@test h2norm(G1-G1) < 1e-6
@test hinfnorm2(G1-G1)[1] < 1e-6
@test nugap(G1, G1)[1] < 1e-6
@test hankelnorm(G1-G1)[1] < 1e-6
@test hinfnorm2(baltrunc2(G1, n=4)[1]-G1)[1] < 1e-10



@test count(1:100) do _
    G1 = ssrand(1,1,3, proper=true)
    G2 = ssrand(1,1,3, proper=true)
    try
        νgap(G1\G2, ss(inv(tf(G1))*tf(G2)))[1] < 1e-8
    catch
        false
    end
end >= 94


##
@testset "unstable baltrunc" begin
    @info "Testing unstable baltrunc"
    for proper = [true, false]
        sys = ssrand(2,3,40; stable=true, proper)
        sysus = ssrand(2,3,2; stable=true, proper)
        sysus.A .*= -1
        sys = sys + sysus

        sysr, hs = RobustAndOptimalControl.baltrunc_coprime(sys, n=20, factorization = RobustAndOptimalControl.DescriptorSystems.glcf)

        @test sysr.nx <= 20
        @test linfnorm(sysr - sys)[1] < 5e-2

        e = poles(sysr)
        @test count(e->real(e)>0, e) == 2 # test that the two unstable poles were preserved
        # bodeplot([sys, sysr])
    end

    ## stab_unstab
    sys = ssrand(2,3,40, stable=false)
    stab, unstab = stab_unstab(sys)
    @test all(real(poles(stab)) .< 0)
    @test all(real(poles(unstab)) .>= 0)
    @test linfnorm2(stab + unstab - sys)[1] < 1e-7

    ## baltrunc_unstab
    sys = ssrand(2,3,40, stable=true)
    sysus = ssrand(2,3,2, stable=true)
    sysus.A .*= -1
    sys = sys + sysus
    sysr, hs = baltrunc_unstab(sys, n=20)
    @test sysr.nx <= 20
    @test linfnorm(sysr - sys)[1] < 1e-2
    # bodeplot([sys, sysr])

end

using RobustAndOptimalControl: baltrunc2, baltrunc_coprime, baltrunc_unstab, _kwarg_names, _select_kwargs
using RobustAndOptimalControl.DescriptorSystems: gsdec, gbalmr, glcf

@testset "keyword-argument forwarding" begin
    @info "Testing keyword-argument forwarding"

    @test _kwarg_names(gsdec) == Set([:prescale, :smarg, :fast, :atol, :atol1, :atol2, :rtol])
    @test :smarg ∈ _kwarg_names(glcf) # glcf slurps its keyword arguments into grcf
    @test :ord ∉ _kwarg_names(gbalmr) # determined by the argument n
    @test _select_kwargs(gsdec, (; smarg=-0.1, atolhsv=1e-3, atol=1e-9)) === (; smarg=-0.1, atol=1e-9)
    @test _select_kwargs(gbalmr, (; smarg=-0.1, atolhsv=1e-3, atol=1e-9)) === (; atolhsv=1e-3, atol=1e-9)

    sys = ssrand(2,3,20, stable=true)
    sysus = ssrand(2,3,2, stable=true)
    sysus.A .*= -1
    sysu = sys + sysus

    # smarg is accepted by gsdec only, atolhsv by gbalmr only
    sysr, _ = baltrunc_unstab(sysu; n=12, smarg=-0.01, atolhsv=1e-8)
    @test sysr.nx <= 12
    @test_throws ArgumentError baltrunc_unstab(sysu; n=12, smrag=-0.01)

    # smarg is accepted by glcf but not by the default factorization gnlcf
    sysr, _ = baltrunc_coprime(sysu; n=12, factorization=glcf, smarg=-0.01, atolhsv=1e-8)
    @test sysr.nx <= 12
    @test_throws ArgumentError baltrunc_coprime(sysu; n=12, smarg=-0.01)
    @test_throws ArgumentError baltrunc2(sys; n=12, smarg=-0.01)

    # smarg reaches the additive decomposition: poles to the right of smarg are treated as unstable and are preserved
    _, unstab = gsdec(dss(sysu); job="stable", smarg=-0.5)
    nx_unstab = size(unstab.A, 1)
    @test nx_unstab == count(p -> real(p) > -0.5, poles(sysu))
    sysr, _ = baltrunc_unstab(sysu; n=nx_unstab, smarg=-0.5)
    @test sysr.nx == nx_unstab
end
