"""Configuration structs"""

@with_kw struct SolverConfig
    solver::supertype(Vern9) = Vern9()
    abstol::Float64 = 0.5e-8
    reltol::Float64 = 5.0e-8
end

@with_kw struct UnitsConfig{M<:Unitful.Unitlike, L<:Unitful.Unitlike, V<:Unitful.Unitlike,
                            A<:Unitful.Unitlike, T<:Unitful.Unitlike, Q<:Unitful.Time,
                            P<:Unitful.Unitlike}
    m::M = u"Msun"
    l::L = u"kpc"
    v::V = u"km/s"
    a::A = u"km/s/Myr"
    τ::T = u"Gyr"
    t::Q = uconvert(τ, 1.0*l/v )
    p::P = v^2
end
