using GalacticDynamics
using Test
using PythonCall

pyimport("sys")."path".append("")
accelerations_py = pyimport("accelerations")

@testset "AllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = GalacticDynamics.AllenSantillanHalo(m, a, Λ, γ)
    pot_py = accelerations_py.AllenSantillan(ustrip.([m, a, Λ, γ])...)
    for i in range(1,100)
        x = 50*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-6
    end
end

@testset "MiyamotoNagaiDisk" begin
    m_gal = 2.325e7*u"Msun"
    m =500.0*m_gal  # Msun
    a = 12.0*u"kpc"     # kpc
    b = 4.0*u"kpc"    # kpc
    pot = GalacticDynamics.MiyamotoNagaiDisk(m, a, b)
    pot_py = accelerations_py.MiyamotoNagai(ustrip.([m, a, b])...)
    for i in range(1,100)
        x = 50*rand(3)*u"kpc"
        @test ustrip.(acceleration(pot,x)) ≈ pyconvert(Vector{Float64},pot_py.accel(ustrip.(x)...)) rtol=5.e-6
    end
end