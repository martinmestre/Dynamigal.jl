@testset "AccelerationsAllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    pot_py = accelerations_py.AllenSantillan(adimensional(m, a, Λ, γ)...)
    for i in range(1,2)
        x = 50*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, ustrip.(x)) rtol=5.e-10
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-6
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
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, ustrip.(x)) rtol=5.e-10
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-6
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
        @test ustrip.(acceleration(pot,x)) ≈ acceleration(pot, ustrip.(x)) rtol=5.e-10
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-6

    end
end

@testset "AccelerationsKepler" begin
    m = 1000𝕦.m  # Msun
    pot = Kepler(m)
    Kepler_accel(pot::Kepler, x::Vector{<:Real}) = -G*pot.m/sqrt(x'x)^3 .* x
    for i in range(1,200)
        x = 50*rand(3)
        @test acceleration(pot, x) ≈ Kepler_accel(pot, x) rtol=5.e-10
    end
end

@testset "ConcentrationNFW" begin
    for i in range(1,200)
        m = rand()*10^12*𝕦.m  # Msun
        a = 20*rand()*𝕦.l
        pot = NFW(m, a)
        c = concentration(pot)
        pot₂ = NFW(m, c)
        @test pot₂.a ≈ pot.a rtol=5.e-10
     end
end