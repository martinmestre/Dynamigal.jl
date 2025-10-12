@testset "ODEsLargeCloudMW" begin
    m_p = 1.0e8
    a_p = 5.0
    m_d = 1.0e10
    a_d = 7.0
    b_d = 2.0
    m_h = 5.0e11
    a_h = 50.0
    mp_array = Vector{MacroParticle}(undef, 2)
    pot₁ = Plummer(m_p, a_p)
    pot₂ = MiyamotoNagaiDisk(m_d, a_d, b_d)
    event₁ = Event(30ones(3), 200ones(3))
    event₂ = Event(-30ones(3), 100ones(3))
    mp_array[1] = MacroParticle(pot₁ + pot₂, event₁)
    mp_array[2] = MacroParticle(pot₁ + pot₂, event₂)
    mps = MacroParticleSystem(mp_array...)
    cloudMW = LargeCloudMW(mps)
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
    u = SA[x...,v...].*rand(12)
    ode₁ = ode(u, mps, 0.0)
    ode₂ = ode(u, cloudMW, 0.0)
    @test ode₁ ≈ ode₂  rtol=5.e-16
    @show ode₁ ode₂
    a = @benchmark ode($u, $mps, 0.0) samples=100 seconds=50
    b = @benchmark ode($u, $cloudMW, 0.0) samples=100 seconds=50
    display(a)
    display(b)
end
