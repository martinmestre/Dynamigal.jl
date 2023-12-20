"""Configuration/option structs and conversion functions"""

"""Solver algorithm"""
@with_kw struct SolverOptions
    abstol::Float64 = 0.5e-10
    reltol::Float64 = 5.0e-10
end
ntSolverOptions() = ntfromstruct(SolverOptions())

@with_kw struct SolverConfig
    ode::supertype(Vern9) = Vern9()
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
code_units(::Nothing) = nothing
code_units(x::L) where {L<:Unitful.Length} = ustrip(uconvert(𝕦.l, x))
code_units(v::V) where {V<:Unitful.Velocity} = ustrip(uconvert(𝕦.v, v))
code_units(t::T) where {T<:Unitful.Time} = uconvert(𝕦.τ, t)/𝕦.t
code_units(x::Vector{L}) where {L<:Unitful.Length} = code_units.(x)
code_units(v::Vector{V}) where {V<:Unitful.Velocity} = code_units.(v)
code_units(x::Vector{L}, t::T) where {L<:Unitful.Length, T<:Union{Unitful.Time,Nothing}} =
    code_units(x), code_units(t)
code_units(x::Vector{L}, v::Vector{V}) where {L<:Unitful.Length, V<:Unitful.Velocity} =
    code_units(x), code_units(v)
code_units(x::Vector{L}, v::Vector{V}, t::T) where {L<:Unitful.Length, V<:Unitful.Velocity, T<:Unitful.Time} =
    code_units(x), code_units(v), code_units(t)


"""Physical units"""
function physical_units(x::T, s::Symbol) where {T<:Real}
    if s==:l
        return x*𝕦.l
    elseif s==:v
        return  x*𝕦.v
    elseif s==:t
        return x*𝕦.t
    elseif s==:a
        return x*𝕦.a
    end
end

physical_units(x::Vector{L}, t::T) where {L<:Real, T<:Real} =
    physical_units.(x,:l), physical_units(t,:t)
physical_units(x::Vector{L}, v::Vector{V}) where {L<:Real, V<:Real} =
    physical_units.(x,:l), physical_units_v.(v,:v)
physical_units(x::Vector{L}, v::Vector{V}, t::T) where {L<:Real, V<:Real,T<:Real} =
    physical_units.(x,:l), physical_units.(v,:v), physical_units(t,:t)


"""Cosmos configuration"""
@with_kw struct CosmosConfig{F<:Unitful.Frequency, D<:Unitful.Density}
    H₀::F = 67.66u"km/Mpc/s"
    ρ_c::D = uconvert(u"Msun/kpc^3", 3H₀^2/(8π*u"G") )
end