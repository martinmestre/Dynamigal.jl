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
        @test sol.x[end] ≈ sol₂.x[end] rtol=5.0e-7
        @test sol.v[end] ≈ sol₂.v[end] rtol=5.0e-7
    end
end

@testset "CircularOrbitsKepler" begin
    m =1𝕦.m  # Msun
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
            @test sol.x[1:2, j] ≈ sol₂[j] rtol=5.0e-7
        end
    end
end



"""Precision test between Gala and GalacticDynamics for the NFW
Adding `atol=0.5e-16` in both Gala and GalacticDynamics improves the precision by two orders of magnitude, only when `rtol` is already very small (`< 10^{-16}`)."""


@testset "SingleOrbitNFWvsGala" begin
    Δt = 0.01
    n_step = 1000
    t₁ = 0.0
    t₂ = t₁ + n_step*Δt
    @show t₂
    t_range = (t₁, t₂)
    x₀ = -50.0*SA[1,0,0]
    v₀ = 200.0*SA[0,1,0]
    m = 10.0^12*𝕦.m  # Msun
    a = 20.0*𝕦.l
    pot = NFW(m, a)
    c = concentration(pot)
    f(x) = log(1+x)-x/(1+x)
    m_g = m/f(c)
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.NFWPotential(Py(adimensional(m_g))*au.Msun, Py(adimensional(a))*au.kpc, units=usys)
    w₀ = gd.PhaseSpacePosition(pos=Py(x₀)*au.kpc, vel=Py(v₀)*au.kpc/au.Gyr)

    orb₁ = pot_Gala.integrate_orbit(w₀, dt=Δt*au.Gyr, t1=t₁, t2=t₂*au.Gyr,
                                Integrator=gi.DOPRI853Integrator,
                                Integrator_kwargs=Py(Dict("rtol"=>5.0e-8)))
    orb₂ = pot_Gala.integrate_orbit(w₀, dt=Δt*au.Gyr, t1=t₁, t2=(t₂)*au.Gyr,
                                    Integrator=gi.DOPRI853Integrator,
                                    Integrator_kwargs=Py(Dict("rtol"=>5.0e-11)))
    orb₃ = pot_Gala.integrate_orbit(w₀, dt=Δt*au.Gyr, t1=t₁, t2=(t₂)*au.Gyr,
                                    Integrator=gi.DOPRI853Integrator,
                                    Integrator_kwargs=Py(Dict("rtol"=>5.0e-20, "atol"=>0.5e-20)))
    orb₄ = evolve(pot, x₀, v₀, t_range, Vern7(); options=ntSolverOptions(;reltol=5.0e-8, saveat=Δt))
    orb₅ = evolve(pot, x₀, v₀, t_range, Vern7(); options=ntSolverOptions(;reltol=5.0e-11,saveat=Δt))
    orb₆ = evolve(pot, x₀, v₀, t_range, Vern7(); options=ntSolverOptions(;reltol=5.0e-16, abstol=0.5e-16, saveat=Δt))
    orb₀ = evolve(pot, x₀, v₀, t_range, Vern7(); options=ntSolverOptions(; saveat=Δt))
    @test orb₄.x[1,end] ≈ pyconvert(Float64,orb₁.x[-1].value)  rtol=5.0e-8
    @test orb₅.x[1,end] ≈ pyconvert(Float64,orb₂.x[-1].value)  rtol=5.0e-10
    @test orb₆.x[1,end] ≈ pyconvert(Float64,orb₃.x[-1].value)  rtol=5.0e-12
    @test orb₀.x[1,end] ≈ pyconvert(Float64,orb₃.x[-1].value)  rtol=5.0e-12
end


@testset "OrbitsNFWvsGala" begin
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    Δt = 0.01
    n_step = 1000
    t₁ = 0.0
    t₂ = t₁ + n_step*Δt
    t_range = (t₁, t₂)
    x₀ = 30*[1,0,1]
    v₀ = 200*[0,1,0]
    m = 10^12*𝕦.m  # Msun
    a = 20*𝕦.l
    pot = NFW(m, a)
    c = concentration(pot)
    f(x) = log(1+x)-x/(1+x)
    m_g = m/f(c)
    pot_Gala = gp.NFWPotential(Py(adimensional(m_g))*au.Msun, Py(adimensional(a))*au.kpc, units=usys)
    @show pot_Gala

    # Gala.py solution
    w₀ = gd.PhaseSpacePosition(pos=Py(x₀)*au.kpc, vel=Py(v₀)*au.kpc/au.Gyr)
    orb_gala = pot_Gala.integrate_orbit(w₀, dt=Δt*au.Gyr, t1=t₁, t2=t₂*au.Gyr,
            Integrator=gi.DOPRI853Integrator,Integrator_kwargs=Py(Dict("rtol"=>5.0e-14, "atol"=>0.5e-14)))
    orb_gala_t = pyconvert(Vector{Float64}, orb_gala.t)
    orb_gala_x = pyconvert(Vector{Float64}, orb_gala.x)
    orb_gala_y = pyconvert(Vector{Float64}, orb_gala.y)
    orb_gala_z = pyconvert(Vector{Float64}, orb_gala.z)
    # GalacticDynamics.jl solution
    sol = evolve(pot, x₀, v₀, t_range, Tsit5(); options=ntSolverOptions(; reltol=5.0e-14, abstol=0.5e-14,saveat=Δt))
    orb_t = ustrip.(physical_units.(sol.t,:t))
    orb_x = sol.x[1,:]

    @test orb_t[end] ≈ orb_gala_t[end] rtol=5.0e-12
    @test orb_x[end] ≈ orb_gala_x[end] rtol=5.0e-12

end

# @testset "ODEAllenSantillanHalo" begin
#     m_gal = 2.325e7*u"Msun"
#     m =1018.0*m_gal  # Msun
#     a = 2.562*u"kpc"     # kpc
#     Λ = 200.0*u"kpc"    # kpc
#     γ = 2.0
#     pot = AllenSantillanHalo(m, a, Λ, γ)
#     for i in range(1,20)
#         w₀ = 50*rand(6)
#         x₀ = w₀[1:3]u"kpc"
#         v₀ = w₀[4:6]u"km/s"
#         t_range = (0.0,100.0).*𝕦.t
#         sol = evolve(pot, x₀, v₀, t_range)
#         sol₂ = evolve(pot, x₀, v₀, t_range; options=ntSolverOptions(; reltol=5.0e-12))
#         @test sol.x[end] ≈ sol₂.x[end] rtol=5.0e-7
#         @test sol.v[end] ≈ sol₂.v[end] rtol=5.0e-7
#     end
# end

@testset "OrbitsMacroParticleSystem" begin
    m = 500
    n = 3
    m_p = 1.0e8
    a_p = 5.0
    m_d = 1.0e10
    a_d = 7.0
    b_d = 2.0
    m_h = 5.0e11
    a_h = 50.0
    pot = Vector{CompositePotential}(undef, n+1)
    mp_array = Vector{MacroParticle}(undef, n+1)
    for i in eachindex(pot)
        pot₁ = Plummer(m_p*rand(), a_p*rand())
        pot₂ = MiyamotoNagaiDisk(m_d*rand(), a_d*rand(), b_d*rand())
        pot₃ = Hernquist(m_h*rand(), a_h*rand())
        pot[i] = CompositePotential(pot₁,pot₂,pot₃)
        event = Event(30rand(3), 200rand(3))
        mp_array[i] = MacroParticle(pot[i], event)
    end
    pot[n+1] =  CompositePotential(Kepler(1.0e7), Kepler(0.1))
    mp_array[n+1] = MacroParticle(pot[n+1])
    mps = MacroParticleSystem(mp_array...)
    t_span = (0., 7.)
    @set_system_trait mps GenSys
    @time orbits = evolve(mps, t_span; options=ntSolverOptions(saveat=0.1))
    @set_system_trait mps GenPerfSys
    @time orbits₂ = evolve(mps, t_span; options=ntSolverOptions(saveat=0.1))
    for i in eachindex(pot)
        @test orbits[i].x ≈ orbits₂[i].x  rtol=1.e-8
        @test orbits[i].v ≈ orbits₂[i].v  rtol=1.e-8
    end
end
