"""Examples"""

"""Plummer example"""
function example_Plummer()
    pot = Plummer(10.0^11*𝕦.m, 10.0*𝕦.l)
    x₀ = [10.0, 0.0, 0.0]*𝕦.l
    v₀ = [0.0,50.0,0.0]*𝕦.v
    t_range = (0.0,10.0).*𝕦.τ
    sol = evolve(pot, x₀, v₀, t_range; options=SolverConfig(reltol=5.0e-12))
    return sol
end

function example_MiyamotoNagai()
    m_gal = 2.325e7*𝕦.m
    m =2856.0*m_gal  # Msun
    a = 4.22*𝕦.l    # kpc
    b =0.292*𝕦.l    # kpc
    pot = MiyamotoNagaiDisk(m, a, b)
    x₀ = [10.0, 0.0, 0.0]*𝕦.l
    v₀ = [0.0,10.0,0.0]*𝕦.v
    t_range = (0.0,10.0).*𝕦.τ
    sol = evolve(pot, x₀, v₀, options=SolverConfig(reltol=5.0e-10))
    return sol
end

function example_AllenSantillan()
    m_gal = 2.325e7*𝕦.m
    m =1018.0*m_gal  # Msun
    a = 2.562*𝕦.l     # kpc
    Λ = 200.0*𝕦.l    # kpc
    γ = 2.0
    pot = AllenSantillanHalo(m, a, Λ, γ)
    w₀ = [3.59746558, 8.24013064, -9.17984456, -58.75537855, -147.5572843, 173.06078831]
    x₀ = [10.0,-8.0,7.0]*𝕦.l
    v₀ = w₀[4:6]*𝕦.v
    t_range = (0.0,10.0).*𝕦.τ
    sol = evolve(pot, x₀, v₀, t_range)
    return sol
end

function example_sum_of_potentials()
    m_gal = 2.325e7*𝕦.m
    m =2856.0*m_gal  # Msun
    a = 4.22*𝕦.l     # kpc
    b =0.292*𝕦.l    # kpc
    pot_mn = MiyamotoNagaiDisk(m, a, b)
    pot_pl = Plummer(10.0^11*𝕦.m, 10.0𝕦.l)
    x₀ = [10.0, 0.0, 0.0]*𝕦.l
    v₀ = [0.0,10.0,0.0]*𝕦.v
    t_range = (0.0,10.0).*𝕦.τ
    sol = evolve(pot_mn+pot_pl, x₀, v₀, t_range)
    return sol
end