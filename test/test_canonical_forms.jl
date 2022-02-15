using RobustAndOptimalControl, ControlSystems
import RobustAndOptimalControl: cdf2rdf, blockdiagonalize

function isblockdiagonal(A)
    complex_inds = findall(diag(A, -1) .!= 0)
    diag(A, -1) ≈ -diag(A, 1) || (return false)
    A = A - diagm(diag(A)) # remove main diagonal
    A = A - diagm(1=>diag(A, 1))
    A = A - diagm(-1=>diag(A, -1))
    all(iszero, A)
end


@testset "modal form" begin
    @info "Testing modal form"

    X = randn(5,5) # odd number to enforce at least one real eigval
    E = eigen(X)
    D,V = E
    Db, Vb = cdf2rdf(E)
    @test isblockdiagonal(Db)
    @test Vb*Db ≈ X*Vb
    @test sort(real(D)) ≈ sort(diag(Db)) # real values on diagonal
    ivals = [diag(Db, -1); diag(Db, 1)] # imag values on 1/-1 diagonals
    @test all(v ∈ ivals for v in imag(D)) 

    Xb, T = blockdiagonalize(X)
    @test T\X*T ≈ Xb
    @test isblockdiagonal(Xb)

    sys = ssrand(1,1,5) # odd number to enforce at least one real eigval
    sysm, T = modal_form(sys)
    @test tf(sys) ≈ tf(sysm)
    @test sysm ≈ similarity_transform(sys, T)

    complex_inds = findall(diag(sysm.A, -1) .!= 0)
    @test all(sysm.C[i] > 0 for i ∈ complex_inds) # test coefficient convention


    # test balanced property
    sys = ssrand(1,1,1, proper=true)
    sysm, T = modal_form(sys)
    @test tf(sys) ≈ tf(sysm)
    @test sysm ≈ similarity_transform(sys, T)
    @test abs(sysm.C[1]) ≈ abs(sysm.B[1])

    # test C1 property
    sys = ssrand(1,1,1, proper=true)
    sysm, T = modal_form(sys, C1=true)
    @test tf(sys) ≈ tf(sysm)
    @test sysm ≈ similarity_transform(sys, T)
    @test sysm.C[1] ≈ 1


    sys = ssrand(2,3,5) # test that it works for MIMO
    sysm, T = modal_form(sys)
    @test tf(sys) ≈ tf(sysm)
    @test isblockdiagonal(sysm.A)

    @test sysm ≈ similarity_transform(sys, T)


    # test with repeated eigenvalues
    X = ss(tf(1,[1,1,1])^2).A
    E = eigen(X)
    Db, Vb = cdf2rdf(E)
    @test isblockdiagonal(Db)
    @test Vb*Db ≈ X*Vb
    @test sort(real(E.values)) ≈ sort(diag(Db)) # real values on diagonal
    ivals = [diag(Db, -1); diag(Db, 1)] # imag values on 1/-1 diagonals
    @test all(v ∈ ivals for v in imag(E.values)) 
end




@testset "schur form" begin
    @info "Testing schur form"

    sys = ssrand(1,1,5) # odd number to enforce at least one real eigval
    sysm, T = schur_form(sys)
    @test tf(sys) ≈ tf(sysm)

    @test sysm ≈ similarity_transform(sys, T)


    sys = ssrand(2,3,5) # test that it works for MIMO
    sysm, _ = schur_form(sys)
    @test tf(sys) ≈ tf(sysm)
    @test all(d->all(iszero, sysm.A[diagind(sysm.A, d)]), -5:-2)
end



@testset "hess form" begin
    @info "Testing hess form"

    sys = ssrand(1,1,5) # odd number to enforce at least one real eigval
    sysm, T = hess_form(sys)
    @test tf(sys) ≈ tf(sysm)

    @test sysm ≈ similarity_transform(sys, T)


    sys = ssrand(2,3,5) # test that it works for MIMO
    sysm, _ = hess_form(sys)
    @test tf(sys) ≈ tf(sysm)
    @test all(d->all(iszero, sysm.A[diagind(sysm.A, d)]), -5:-2)
end


