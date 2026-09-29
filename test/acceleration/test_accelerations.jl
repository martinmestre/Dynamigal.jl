@testset "AccelerationsAllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    pot_py = accelerations_py.AllenSantillan(adimensional(m, a, Λ)...,γ)
    for i in range(1,20)
        x = 200*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, adimensional(x)) rtol=5.e-14
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(adimensional(x)...)) rtol=5.e-14
    end
end

@testset "AccelerationsMiyamotoNagai" begin
    m_gal = 2.325e7*u"Msun"
    m =500.0*m_gal  # Msun
    a = 12.0*u"kpc"     # kpc
    b = 4.0*u"kpc"    # kpc
    pot = MiyamotoNagai(m, a, b)
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
    for i in range(1,2)
        m = rand()*10^12*𝕦.m  # Msun
        a = 20*rand()*𝕦.l
        pot = NFW(m, a)
        pot₂ = NFW(:m, :c, m, pot.c)
        @show pot pot₂
        @test pot₂.a ≈ pot.a rtol=5.e-14
     end
end

@testset "AccelerationMilkyWayBovy2014vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.BovyMWPotential2014()
    pot = MilkyWayBovy2014()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        acc_gala = pot_Gala.acceleration(x).to(au.kpc/au.Gyr^2).value
        for j in 1:3
            @test acceleration(pot,x)[j] ≈ pyconvert(Float64,acc_gala[j-1][0]) rtol=5.0e-11
        end
    end
end
@testset "AccelerationMilkyWayPriceWhelan2017vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.MilkyWayPotential(version="v1")
    pot = MilkyWayPriceWhelan2017()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        acc_gala = pot_Gala.acceleration(x).to(au.kpc/au.Gyr^2).value
        for j in 1:3
            @test acceleration(pot,x)[j] ≈ pyconvert(Float64,acc_gala[j-1][0]) rtol=5.0e-11
        end
    end
end

@testset "AccelerationMilkyWayPriceWhelan2022vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.MilkyWayPotential(version="v2")
    pot = MilkyWayPriceWhelan2022()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        acc_gala = pot_Gala.acceleration(x).to(au.kpc/au.Gyr^2).value
        for j in 1:3
            @test acceleration(pot,x)[j] ≈ pyconvert(Float64,acc_gala[j-1][0]) rtol=5.0e-11
        end
    end
end

@testset "AccelerationsMacroParticleSystem" begin
    m = 5
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
        pot₂ = MiyamotoNagai(m_d*rand(), a_d*rand(), b_d*rand())
        pot₃ = Hernquist(m_h*rand(), a_h*rand())
        pot[i] = CompositePotential(pot₁,pot₂,pot₃)
        event = Event(30rand(3), 200rand(3))
        mp_array[i] = MacroParticle(pot[i], event)
    end
    pot[n+1] =  CompositePotential(Kepler(1.0e7), Kepler(0.1))
    mp_array[n+1] = MacroParticle(pot[n+1])
    mps = MacroParticleSystem(mp_array...)
    # @show typeof(mps) typeof(mps.macroparticles) typeof(mps.accelerations)
    x = vcat([[mps[i].event.x for i ∈ eachindex(mps)]...]...)
    v = vcat([[mps[i].event.v for i ∈ eachindex(mps)]...]...)
    for j = 1:2
        u = SA[x...,v...].*rand(6*(n+1))
        acc = acceleration!(mps,u)
        acc₂ = acceleration_c!(mps,u)
        @test acc ≈ acc₂ rtol=5.e-14
        a = @benchmark acceleration!($mps,$u) samples=100 seconds=1000
        b = @benchmark acceleration_c!($mps,$u) samples=100 seconds=1000
        println("j = $j")
        display(a)
        display(b)
    end
    for Δx ∈ [0.01, 0.001, 0.0001, 0.00001, 0.000001]
        x = Δx*[1,0,0]
        y = zeros(MVector{3*(n+1), Float64})  # I am not using the true mps positions.
        # @show acceleration(mps, y, x)
        c = @benchmark acceleration($mps, $y, $x) samples=100 seconds=1000
        println("Δx = $(Δx)")
        display(c)
    end
end

@testset "AccelerationsLargeCloudMW" begin
    m_p = 1.0e8
    a_p = 5.0
    m_d = 1.0e10
    a_d = 7.0
    b_d = 2.0
    m_h = 5.0e11
    a_h = 50.0
    mp_array = Vector{MacroParticle}(undef, 2)
    pot₁ = Plummer(m_p, a_p)
    pot₂ = MiyamotoNagai(m_d, a_d, b_d)
    event₁ = Event(30ones(3), 200ones(3))
    event₂ = Event(-30ones(3), 100ones(3))
    mp_array[1] = MacroParticle(pot₁ + pot₂, event₁)
    mp_array[2] = MacroParticle(pot₁ + pot₂, event₂)
    mps = MacroParticleSystem(mp_array...)
    cloudMW = LargeCloudMW(mps)
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
    u = SA[x...,v...]
    fric_galpy = GalpyFriction(mₚ=1.e-1, rₚ=5.0, σₕ=r->150.0)
    fric_agama = AgamaFriction(mₚ=1.e-1, σₕ=r->150.0)
    fric_tango = TangoFriction(mₚ=1.e-1, σₕ=r->150.0)
    @show fric_galpy fric_agama fric_tango
    acc_c = acceleration_c!(mps,x)
    acc = acceleration!(mps,x)
    acc₂ = acceleration(cloudMW,u)
    acc₃ = acceleration(fric_galpy, cloudMW, u)
    acc₄ = acceleration(fric_agama, cloudMW, u)
    acc₅ = acceleration(fric_tango, cloudMW, u)
    @test acc ≈ acc_c rtol=5.e-14
    @test acc ≈ acc₂ rtol=5.e-14
    @test acc₂ ≈ acc₃ rtol=5.e-14
    @test acc₃ ≈ acc₄ rtol=5.e-14
    @test acc₄ ≈ acc₅ rtol=5.e-14
    @show acc acc_c acc₂ acc₃  acc₄ acc₅
    a_c = @benchmark acceleration_c!($mps,$x) samples=100 seconds=50
    a = @benchmark acceleration!($mps,$x) samples=100 seconds=50
    b = @benchmark acceleration($cloudMW,$u) samples=100 seconds=50
    c = @benchmark acceleration($fric_galpy, $cloudMW,$u) samples=100 seconds=50
    d = @benchmark acceleration($fric_agama, $cloudMW,$u) samples=100 seconds=50
    e = @benchmark acceleration($fric_tango, $cloudMW,$u) samples=100 seconds=50
    display(a_c)
    display(a)
    display(b)
    display(c)
end

@testset "DragsLargeCloudMW" begin
    m_p = 1.0e8
    a_p = 5.0
    m_d = 1.0e10
    a_d = 7.0
    b_d = 2.0
    m_h = 5.0e11
    a_h = 50.0
    mp_array = Vector{MacroParticle}(undef, 2)
    pot₁ = Plummer(m_p, a_p)
    pot₂ = MiyamotoNagai(m_d, a_d, b_d)
    event₁ = Event(30ones(3), 200ones(3))
    event₂ = Event(-30ones(3), 100ones(3))
    mp_array[1] = MacroParticle(pot₁ + pot₂, event₁)
    mp_array[2] = MacroParticle(pot₁ + pot₂, event₂)
    mps = MacroParticleSystem(mp_array...)
    cloudMW = LargeCloudMW(mps)
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
    u = SA[x...,v...]
    fric_galpy = GalpyFriction(mₚ=1.e5, rₚ=5.0, σₕ=r->150.0)
    fric_agama = AgamaFriction(mₚ=1.e5, σₕ=r->150.0)
    fric_tango = TangoFriction(mₚ=1.e5, σₕ=r->150.0)
    @show fric_galpy fric_agama fric_tango
    acc_g = acceleration(fric_galpy, cloudMW, u)
    acc_a = acceleration(fric_agama, cloudMW, u)
    acc_t = acceleration(fric_tango, cloudMW, u)
    @test acc_g ≈ acc_a rtol=5.e-10
    @test acc_g ≈ acc_t rtol=5.e-10
    @test acc_a ≈ acc_t rtol=5.e-10

    @show acc_g acc_a acc_t
    a = @benchmark acceleration($fric_galpy, $cloudMW,$u) samples=100 seconds=50
    b = @benchmark acceleration($fric_galpy, $cloudMW,$u) samples=100 seconds=50
    c = @benchmark acceleration($fric_galpy, $cloudMW,$u) samples=100 seconds=50
    display(a)
    display(b)
    display(c)
end