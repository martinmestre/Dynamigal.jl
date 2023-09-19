"""Examples"""

"""Plummer example"""
function example_Plummer()
    plum = Plummer(10.0^11*u"Msun",10.0u"kpc")
    x₀ = [10.0, 0.0, 0.0]u"kpc"
    v₀ = [0.0,50.0,0.0]u"km/s"
    t_range = (0.0,10.0).*u_T
    @show t_range
    sol = evolve(plum, x₀, v₀, t_range)
    return sol
end