@testset "MilkyWayBovy2014vsGala" begin
    n = 10
    usys = gu.UnitSystem(au.kpc, au.Gyr, au.Msun, au.radian, au.kpc/au.Gyr, au.kpc/au.Gyr^2)
    pot_Gala = gp.BovyMWPotential2014()
    pot = MilkyWayBovy2014()
    for i in 1:n
        s = 0.5*(2rand(3).-1)
        x = 100*s
        @test density(pot,x) ≈ pyconvert(Float64,pot_Gala.density(x).value[0]) rtol=5.0e-11
        acc_gala = pot_Gala.acceleration(x).to(au.kpc/au.Gyr^2).value
        for j in 1:3
            @test acceleration(pot,x)[j] ≈ pyconvert(Float64,acc_gala[j-1][0]) rtol=5.0e-11
        end
    end
end
