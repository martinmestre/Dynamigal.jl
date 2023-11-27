@testset "OrbitsAllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    for i in range(1,20)
        w₀ = 50*rand(6)
        x₀ = w₀[1:3]u"kpc"
        v₀ = w₀[4:6]u"km/s"
        t_range = (0.0,100.0).*𝕦.t
        sol = evolve(pot, x₀, v₀, t_range)
        sol₂ = evolve(pot, x₀, v₀, t_range; options=SolverConfig(reltol=5.0e-12))
        @test sol.x[end] ≈ sol₂.x[end] rtol=5.0e-6
        @test sol.v[end] ≈ sol₂.v[end] rtol=5.0e-6
    end
end

@testset "CircularOrbitsKepler" begin
    m =1.0𝕦.m  # Msun
    pot = Kepler(m)
    function evolve_Kepler_circular(pot, x₀, t)
        a = sqrt(x₀'x₀)
        T = 2π√(a^3/(G*pot.m))
        n =  2π/T
        x = [ a*[cos(n*t[i]), sin(n*t[i])] for i ∈ eachindex(t)]
        return x
    end
    for i in range(1,20)
        r = 10.0*rand()
        x₀ = [r, 0.0, 0.0]
        speed = circular_velocity(pot, x₀)
        v₀ = [0.0, speed, 0.0]
        t_range = (0.0, 10.0)
        sol = evolve(pot, x₀, v₀, t_range)
        sol₂ = evolve_Kepler_circular(pot, x₀, sol.t)
        for j ∈ eachindex(sol.t)
            @test sol.x[1, j] ≈ sol₂[j][1] rtol=5.0e-6
        end
    end
end