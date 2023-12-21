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
        sol₂ = evolve(pot, x₀, v₀, t_range; options=ntSolverOptions(; reltol=5.0e-12))
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
            @test sol.x[1:2, j] ≈ sol₂[j] rtol=5.0e-6
        end
    end
end

@testset "OrbitsNFWvsGala" begin
    usys = gu.UnitSystem(au.kpc, au.Myr, au.Msun, au.radian, au.km/au.s, au.km/au.s/au.Myr)
    t₁, t₂ = 0.0, 10.0
    t_range = code_units.((t₁, t₂).*𝕦.τ)
    Δt = 0.1
    @show t_range
    @show code_units(Δt*𝕦.τ)
    x₀ = 30ones(3)
    v₀ = 100ones(3)
    m = 10^12*𝕦.m  # Msun
    a = 20*𝕦.l
    pot_Gala = gp.NFWPotential(Py(ustrip(m))*au.Msun, Py(ustrip(a))*au.kpc, units=gu.galactic)
    pot = NFW(m, a)
    @show pot_Gala
    for i in range(0,1)
        # Gala.py solution
        w₀ = gd.PhaseSpacePosition(pos=Py(x₀)*au.kpc, vel=Py(v₀)*au.km/au.s)
        orb_gala = pot_Gala.integrate_orbit(w₀, dt=Δt*au.Gyr, t1=t₁, t2=(t₂+Δt)*au.Gyr )
        orb_gala_t = pyconvert(Vector{Float64}, orb_gala.t)
        orb_gala_x = pyconvert(Vector{Float64}, orb_gala.x)
        orb_gala_y = pyconvert(Vector{Float64}, orb_gala.y)
        orb_gala_z = pyconvert(Vector{Float64}, orb_gala.z)
        # GalacticDynamics.jl solution
        sol = evolve(pot, x₀, v₀, t_range; options=ntSolverOptions(; saveat=code_units(Δt*𝕦.τ)))
        @show sol.t
        orb_t = ustrip.(physical_units.(sol.t,:t))*10^3
        orb_x = sol.x[1,:]
        @show orb_t
        @show orb_gala_t
        @show orb_x
        @show orb_gala_x
        @test orb_t[50] ≈ orb_gala_t[50] rtol=5.0e-6
        @test orb_x[50] ≈ orb_gala_x[50] rtol=5.0e-6
    end
end