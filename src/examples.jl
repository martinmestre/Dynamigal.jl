"""Examples"""

"""Plummer example"""
function example_Plummer()
    pot = Plummer(10.0^11*u"Msun",10.0u"kpc")
    x₀ = [10.0, 0.0, 0.0]u"kpc"
    v₀ = [0.0,50.0,0.0]u"km/s"
    t_range = (0.0,10.0).*u_T
    sol = evolve(pot, x₀, v₀, t_range; options=SolverConfig(reltol=5.0e-12))
    @show sol[begin]
    return sol
end

function example_MiyamotoNagai()
    m_gal = 2.325e7*u"Msun"
    m =2856.0*m_gal  # Msun
    a = 4.22*u"kpc"     # kpc
    b =0.292*u"kpc"    # kpc
    pot = MiyamotoNagaiDisk(m, a, b)
    x₀ = [10.0, 0.0, 0.0]u"kpc"
    v₀ = [0.0,10.0,0.0]u"km/s"
    t_range = (0.0,10.0).*u_T
    sol = evolve(pot, x₀, v₀, t_range)
    return sol
end

function example_AllenSantillan()
    m_gal = 2.325e7*u"Msun"
    m =1018.0*m_gal  # Msun
    a = 2.562*u"kpc"     # kpc
    Λ = 200.0*u"kpc"    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    w₀ = [3.59746558, 8.24013064, -9.17984456, -58.75537855, -147.5572843, 173.06078831]
    x₀ = [10.0,-8.0,7.0]u"kpc"
    v₀ = w₀[4:6]u"km/s"
    t_range = (0.0,10.0).*u_T
    sol = evolve(pot, x₀, v₀, t_range)
    return sol
end

function example_sum_of_potentials()
    m_gal = 2.325e7*u"Msun"
    m =2856.0*m_gal  # Msun
    a = 4.22*u"kpc"     # kpc
    b =0.292*u"kpc"    # kpc
    pot_mn = MiyamotoNagaiDisk(m, a, b)
    pot_pl = Plummer(10.0^11*u"Msun",10.0u"kpc")
    x₀ = [10.0, 0.0, 0.0]u"kpc"
    v₀ = [0.0,10.0,0.0]u"km/s"
    t_range = (0.0,10.0).*u_T
    sol = evolve(pot_mn+pot_pl, x₀, v₀, t_range)
    return sol
end