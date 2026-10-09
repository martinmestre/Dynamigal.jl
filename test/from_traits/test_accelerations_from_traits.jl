

@testset "AccelerationsFromPotentialTrait" begin
    function loop(pot, trait, n, δᵣ)
        for i in range(1,n)
            x = 200*rand(3)*u"kpc"
            t = 10.0*u"Gyr"
            r  = sqrt(  dot(x,x)  )

            @test ustrip.(acceleration(trait, pot, x, t)) ≈ ustrip.(acceleration(trait, pot, x)) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, x)) ≈ ustrip.(acceleration(trait, pot, r)*x/r) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, x)) ≈ ustrip.(acceleration(pot, x)) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, r)) ≈ ustrip.(acceleration(pot, r)) rtol=δᵣ

            x, r, t = adimensional(x,r,t)
            @test acceleration(trait, pot, x, t) ≈ acceleration(trait, pot, x) rtol=δᵣ
            @test acceleration(trait, pot, x) ≈ acceleration(trait, pot, r)*x/r rtol=δᵣ
            @test acceleration(trait, pot, x) ≈ acceleration(pot, x) rtol=δᵣ
            @test acceleration(trait, pot, r) ≈ acceleration(pot, r) rtol=δᵣ
        end
    end
    trait = FromPotentialTrait()
    n = 5
    δᵣ = 5.0e-14

    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    loop(pot, trait, n, δᵣ)

    m_gal = 2.325e7*u"Msun"
    m = 1000.0*m_gal  # Msun
    a = 2.0*u"kpc"     # kpc
    pot = Plummer(m, a)
     loop(pot, trait, n, δᵣ)

    # pot = PowerLawCutoff(m=4501365375.06545*u"Msun", α=1.8, c=1.0*u"kpc")
    # loop(pot, trait, n, δᵣ)

    m=4501365375.06545*u"Msun"
    a=40u"kpc"
    pot = NFW(m, a)
    loop(pot, trait, n, δᵣ)

    m=4501365375.06545*u"Msun"
    a=40u"kpc"
    pot = Hernquist(m, a)
    loop(pot, trait, n, δᵣ)

end


@testset "AccelerationsFromPotentialTrait_Aspherical" begin
    function loop(pot, trait, n, δᵣ)
        for i in range(1,n)
            x = 200*rand(3)*u"kpc"
            t = 10.0*u"Gyr"
            r  = sqrt(  dot(x,x)  )

            @test ustrip.(acceleration(trait, pot, x, t)) ≈ ustrip.(acceleration(trait, pot, x)) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, x)) ≈ ustrip.(acceleration(pot, x)) rtol=δᵣ

            x, r, t = adimensional(x,r,t)
            @test acceleration(trait, pot, x, t) ≈ acceleration(trait, pot, x) rtol=δᵣ
            @test acceleration(trait, pot, x) ≈ acceleration(pot, x) rtol=δᵣ
        end
    end
    trait = FromPotentialTrait()
    n = 5
    δᵣ = 5.0e-14


    m=4501365375.06545*u"Msun"
    a=40u"kpc"
    b=3u"kpc"
    pot = MiyamotoNagai(m, a, b)
    loop(pot, trait, n, δᵣ)
end

@testset "AccelerationsFromMassTrait" begin
     function loop(pot, trait, n, δᵣ)
        for i in range(1,n)
            x = 200*rand(3)*u"kpc"
            t = 10.0*u"Gyr"
            r  = sqrt(  dot(x,x)  )

            @test ustrip.(acceleration(trait, pot, x, t)) ≈ ustrip.(acceleration(trait, pot, x)) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, x)) ≈ ustrip.(acceleration(trait, pot, r)*x/r) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, x)) ≈ ustrip.(acceleration(pot, x)) rtol=δᵣ
            @test ustrip.(acceleration(trait, pot, r)) ≈ ustrip.(acceleration(pot, r)) rtol=δᵣ

            x, r, t = adimensional(x,r,t)
            @test acceleration(trait, pot, x, t) ≈ acceleration(trait, pot, x) rtol=δᵣ
            @test acceleration(trait, pot, x) ≈ acceleration(trait, pot, r)*x/r rtol=δᵣ
            @test acceleration(trait, pot, x) ≈ acceleration(pot, x) rtol=δᵣ
            @test acceleration(trait, pot, r) ≈ acceleration(pot, r) rtol=δᵣ
        end
    end
    trait = FromMassTrait()
    n = 5
    δᵣ = 5.0e-14


    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    loop(pot, trait, n, δᵣ)

    m_gal = 2.325e7*u"Msun"
    m = 1000.0*m_gal  # Msun
    a = 2.0*u"kpc"     # kpc
    pot = Plummer(m, a)
    loop(pot, trait, n, δᵣ)

    pot = PowerLawCutoff(m=4501365375.06545*u"Msun", α=1.8, c=1.0*u"kpc")
    loop(pot, trait, n, δᵣ)

    m=4501365375.06545*u"Msun"
    a=40u"kpc"
    pot = NFW(m, a)
    loop(pot, trait, n, δᵣ)

    m=4501365375.06545*u"Msun"
    a=40u"kpc"
    pot = Hernquist(m, a)
    loop(pot, trait, n, δᵣ)
end

