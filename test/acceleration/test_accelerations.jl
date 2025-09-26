@testset "AccelerationsAllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    pot_py = accelerations_py.AllenSantillan(adimensional(m, a, Λ)...,γ)
    for i in range(1,100)
        x = 50*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, adimensional(x)) rtol=5.e-14
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(adimensional(x)...)) rtol=5.e-14
    end
end

@testset "AccelerationsMiyamotoNagaiDisk" begin
    m_gal = 2.325e7*u"Msun"
    m =500.0*m_gal  # Msun
    a = 12.0*u"kpc"     # kpc
    b = 4.0*u"kpc"    # kpc
    pot = MiyamotoNagaiDisk(m, a, b)
    pot_py = accelerations_py.MiyamotoNagai(ustrip.([m, a, b])...)
    for i in range(1,100)
        x = 50*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, ustrip.(x)) rtol=5.e-14
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-14
    end
end

@testset "AccelerationsPlummer" begin
    m_gal = 2.325e7*u"Msun"
    m = 1000.0*m_gal  # Msun
    a = 2.0*u"kpc"     # kpc
    pot = Plummer(m, a)
    pot_py = accelerations_py.Plummer(ustrip.([m, a])...)
    for i in range(1,100)
        x = 50*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, ustrip.(x)) rtol=5.e-14
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-14

    end
end

@testset "AccelerationsKepler" begin
    m = 1000𝕦.m  # Msun
    pot = Kepler(m)
    Kepler_accel(pot::Kepler, x::Vector{<:Real}) = -G*pot.m/sqrt(x'x)^3 .* x
    for i in range(1,200)
        x = 50*rand(3)
        @test acceleration(pot, x) ≈ Kepler_accel(pot, x) rtol=5.e-14
    end
end

@testset "ConcentrationNFW" begin
    for i in range(1,200)
        m = rand()*10^12*𝕦.m  # Msun
        a = 20*rand()*𝕦.l
        pot = NFW(m, a)
        c = concentration(pot)
        pot₂ = NFW(m, c)
        @test pot₂.a ≈ pot.a rtol=5.e-14
     end
end




@testset "AccelerationsMacroParticleSystem" begin
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
    SystemTrait(::Type{typeof(mps)}) = GenPerfSys()
    orbits = evolve(mps, t_span)
    @show SystemTrait(typeof(mps))
    # @set_trait typeof(mps) GenPerfSys()
    # orbits₂ = evolve(mps, t_span)
end
