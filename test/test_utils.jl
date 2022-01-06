G = ssrand(2,3,4)

show_construction(G)

for ny = 1:3, nu = 1:3, nx = 1:3, proper=(true, false)
    G = ssrand(ny, nu, nx; proper)
    v = vec(G)
    G2 = vec2sys(v, ny, nu)
    @test G2 == G
end