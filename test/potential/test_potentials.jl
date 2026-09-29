@testset "PotentialMilkyWayBovy2014vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.BovyMWPotential2014()
    pot = MilkyWayBovy2014()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        @test density(pot,x) ≈ pyconvert(Float64,pot_Gala.density(x).value[0]) rtol=5.0e-11
        conv_Gala = uconvert(𝕦.p, pyconvert(Float64, pot_Gala.energy(x).value[0])*u"kpc^2/Myr^2")
        @test potential(pot,x) ≈ ustrip(conv_Gala) rtol=5.0e-11
    end
end


@testset "PotentialMilkyWayPriceWhelan2017vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.MilkyWayPotential(version="v1")
    pot = MilkyWayPriceWhelan2017()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        @test density(pot,x) ≈ pyconvert(Float64,pot_Gala.density(x).value[0]) rtol=5.0e-11
        conv_Gala = uconvert(𝕦.p, pyconvert(Float64, pot_Gala.energy(x).value[0])*u"kpc^2/Myr^2")
        @test potential(pot,x) ≈ ustrip(conv_Gala) rtol=5.0e-11
    end
end
@testset "PotentialMilkyWayPriceWhelan2022vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.MilkyWayPotential(version="v2")
    pot = MilkyWayPriceWhelan2022()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        @test density(pot,x) ≈ pyconvert(Float64,pot_Gala.density(x).value[0]) rtol=5.0e-11
        conv_Gala = uconvert(𝕦.p, pyconvert(Float64, pot_Gala.energy(x).value[0])*u"kpc^2/Myr^2")
        @test potential(pot,x) ≈ ustrip(conv_Gala) rtol=5.0e-11
    end
end