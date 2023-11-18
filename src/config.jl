"""Configuration structs"""

@with_kw struct SolverConfig
    solver::supertype(Vern9) = Vern9()
    abstol::Float64 = 0.5e-8
    reltol::Float64 = 5.0e-8
end

@with_kw struct UnitsConfig{U,Q} where {U<:Unitful.Unitlike,Q<:Unitful.Time}
    m::U = u"Msun"
    l::U = u"kpc"
    v::U = u"km/s"
    a::U = u"km/s^2"
    t::Q = uconvert(u"Gyr", 1*l/v)
    p::U = v^2
    τ::U = u"Gyr"
end
