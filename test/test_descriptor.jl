using RobustAndOptimalControl, ControlSystems
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
