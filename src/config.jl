"""Configuration structs and conversion functions"""

"""Solver Configuration"""
@with_kw struct SolverConfig
    solver::supertype(Vern9) = Vern9()
    abstol::Float64 = 0.5e-8
    reltol::Float64 = 5.0e-8
end

"""Units configuration"""
@with_kw struct UnitsConfig{M<:Unitful.Unitlike, L<:Unitful.Unitlike, V<:Unitful.Unitlike,
                            A<:Unitful.Unitlike, T<:Unitful.Unitlike, Q<:Unitful.Time,
                            P<:Unitful.Unitlike}
    m::M = u"Msun"
    l::L = u"kpc"
    v::V = u"km/s"
    a::A = u"km/s/Myr"
    τ::T = u"Gyr"
    t::Q = uconvert(τ, 1.0*l/v ) # This is the code natural time unit corresponding to G.
    p::P = v^2
end


"""Code units"""
code_units(x::Vector{L}) where {L<:Unitful.Length} = ustrip(uconvert.(𝕦.l, x))
code_units(v::Vector{V}) where {V<:Unitful.Velocity} = ustrip(uconvert.(𝕦.v, v))
code_units(t::T) where {T<:Unitful.Time} = uconvert(𝕦.τ, t)/𝕦.t
code_units(x::Vector{L}, t::T) where {L<:Unitful.Length, T<:Unitful.Time} =
    code_units(x), code_units(t)
code_units(x::Vector{L}, v::Vector{V}) where {L<:Unitful.Length, V<:Unitful.Velocity} =
    code_units(x), code_units(v)
code_units(x::Vector{L}, v::Vector{V}, t::T) where {L<:Unitful.Length, V<:Unitful.Velocity,T<:Unitful.Time} =
    code_units(x), code_units(v), code_units(t)


"""Physical units"""
physical_units(x::Vector{L}) where {L<:Real} = x*𝕦.l
physical_units(v::Vector{V}) where {V<:Real} = v*𝕦.v
physical_units(t::T) where {T<:Real} = t*𝕦.t
physical_units(x::Vector{L}, t::T) where {L<:Real, T<:Real} =
    physical_units(x), physical_units(t)
physical_units(x::Vector{L}, v::Vector{V}) where {L<:Real, V<:Real} =
    physical_units(x), physical_units(v)
physical_units(x::Vector{L}, v::Vector{V}, t::T) where {L<:Real, V<:Real,T<:Real} =
    physical_units(x), physical_units(v), physical_units(t)