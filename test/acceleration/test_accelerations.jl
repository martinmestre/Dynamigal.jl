@testset "AllenSantillanHalo" begin
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    pot_py = accelerations_py.AllenSantillan(m, a, Λ, γ)



end