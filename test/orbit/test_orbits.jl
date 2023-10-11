@testset "OrbitsAllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = GalacticDynamics.AllenSantillanHalo(m, a, Λ, γ)
    for i in range(1,2)
        w₀ = 50*rand(6)
        x₀ = w₀[1:3]u"kpc"
        v₀ = w₀[4:6]u"km/s"
        t_range = (0.0,100.0).*GalacticDynamics.u_T
        sol = evolve(pot, x₀, v₀, t_range)
        sol₂ = evolve(pot, x₀, v₀, t_range; options=SolverConfig(reltol=5.0e-12))
        @test sol.x[end] ≈ sol₂.x[end] rtol=5.0e-6
        @test sol.v[end] ≈ sol₂.v[end] rtol=5.0e-6
    end
end

