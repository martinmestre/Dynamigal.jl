@testset "Orbits@AllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    for i in range(1,100)
        w₀ = 50*rand(6)*u"kpc"
        x₀ = w₀[1:3]u"kpc"
        v₀ = w₀[4:6]u"km/s"
        t_range = (0.0,10.0).*u_T
        sol = evolve(pot, x₀, v₀, t_range)
        sol₂ = evolve(pot, x₀, v₀, t_range; options=SolverConfig(abstol=0.5e-5, reltol=5.0e-5))
        @test sol ≈ sol₂ reltol=5.0e-6
    end
end

