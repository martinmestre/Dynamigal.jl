@testset "PotentialsFromAccelerationTrait" begin
    function loop(pot, trait, n, δᵣ)
        for i in range(1, n)
            x = 200*rand(3)*u"kpc"
            t = 10.0*u"Gyr"
            r  = sqrt(  dot(x,x)  )

            @test ustrip(potential(trait, pot, x, t)) ≈ ustrip(potential(trait, pot, x)) rtol=δᵣ
            @test ustrip(potential(trait, pot, x)) ≈ ustrip(potential(trait, pot, r)) rtol=δᵣ
            @test ustrip(potential(trait, pot, x)) ≈ ustrip(potential(pot, x)) rtol=δᵣ
            @test ustrip(potential(trait, pot, r)) ≈ ustrip(potential(pot, r)) rtol=δᵣ

            x, r, t = adimensional(x,r,t)
            @test potential(trait, pot, x, t) ≈ potential(trait, pot, x) rtol=δᵣ
            @test potential(trait, pot, x) ≈ potential(trait, pot, r) rtol=δᵣ
            @test potential(trait, pot, x) ≈ potential(pot, x) rtol=δᵣ
            @test potential(trait, pot, r) ≈ potential(pot, r) rtol=δᵣ
        end
    end

    δᵣ = 5.0e-12
    n = 5
    trait = FromAccelerationTrait()

    # Uses AD for acceleration, while I do not write the acceleration function.
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

    # Uses FromAccelerationTrait as default, while I do not write the potential function.
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

@testset "PotentialsFromDensityTrait" begin
    function loop(pot, trait, n, δᵣ)
        for i in range(1, n)
            x = 200*rand(3)*u"kpc"
            t = 10.0*u"Gyr"
            r  = sqrt(  dot(x,x)  )

            @test ustrip(potential(trait, pot, x, t)) ≈ ustrip(potential(trait, pot, x)) rtol=δᵣ
            @test ustrip(potential(trait, pot, x)) ≈ ustrip(potential(trait, pot, r)) rtol=δᵣ
            @test ustrip(potential(trait, pot, x)) ≈ ustrip(potential(pot, x)) rtol=δᵣ
            @test ustrip(potential(trait, pot, r)) ≈ ustrip(potential(pot, r)) rtol=δᵣ

            x, r, t = adimensional(x,r,t)
            @test potential(trait, pot, x, t) ≈ potential(trait, pot, x) rtol=δᵣ
            @test potential(trait, pot, x) ≈ potential(trait, pot, r) rtol=δᵣ
            @test potential(trait, pot, x) ≈ potential(pot, x) rtol=δᵣ
            @test potential(trait, pot, r) ≈ potential(pot, r) rtol=δᵣ
        end
    end

    δᵣ = 5.0e-11
    n = 5
    trait = FromDensityTrait()

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

    # Uses FromAccelerationTrait as default, while I do not write the potential function.
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