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
end >= 95

