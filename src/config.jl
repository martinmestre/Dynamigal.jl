"""Configuration/option structs and conversion functions"""


"""Solver algorithm"""
@with_kw struct SolverOptions  <:AbstractConfig
    abstol::Float64 = 0.5e-12
    reltol::Float64 = 5.0e-12
end
ntSolverOptions(; kwargs...) = (; ntfromstruct(SolverOptions())..., kwargs...)

@with_kw struct SolverConfig  <:AbstractConfig
    ode::supertype(Vern9) = Vern9()
end


"""Units configuration"""
@with_kw struct UnitsConfig{M<:Unitful.Unitlike, L<:Unitful.Unitlike, T<:Unitful.Unitlike,
                            V<:Unitful.Unitlike, A<:Unitful.Unitlike,P<:Unitful.Unitlike,
                            V₁<:Unitful.Unitlike, A₁<:Unitful.Unitlike, R<:Unitful.Unitlike}  <:AbstractConfig
    m::M = u"Msun" # ver si cambio a Giga Msun
    l::L = u"kpc"
    t::T = u"Gyr"
    v::V = l/t    # This is the code unit for velocity
    ν::V₁ = u"km/s"  # Just for IC and display
    a::A = v/t    # This is the code unit for acceleration
    α::A₁ = u"km/s/Myr"  # Just for IC and display
    p::P = v^2  # Code unit for potentials
    ρ::R = m/l^3
end



"""Code units"""
code_units(::Nothing) = nothing
code_units(x::L) where {L<:Unitful.Length} = uconvert(𝕦.l, x)
code_units(v::V) where {V<:Unitful.Velocity} = uconvert(𝕦.v, v)
code_units(t::T) where {T<:Unitful.Time} = uconvert(𝕦.t, t)
code_units(a::A) where {A<:Unitful.Acceleration} = uconvert(𝕦.a, a)
code_units(m::M) where {M<:Unitful.Mass} = uconvert(𝕦.m, m)
code_units(x::Vector{L}) where {L<:Unitful.Length} = code_units.(x)
code_units(v::Vector{V}) where {V<:Unitful.Velocity} = code_units.(v)
code_units(a::Vector{A}) where {A<:Unitful.Acceleration} = code_units.(a)
code_units(m::Vector{M}) where {M<:Unitful.Mass} = code_units.(m)
code_units(x...) = code_units.([x...])


adimensional(x::T) where {T} = ustrip.(code_units(x))
adimensional(x...) = adimensional.([x...])

"""Physical units"""
function physical_units(x::T, s::Symbol) where {T<:Real}
    if s==:l
        return x*𝕦.l
    elseif s==:t
        return x*𝕦.t
    elseif s==:v
        return x*𝕦.v
    elseif s==:a
        return x*𝕦.a
    elseif s==:m
        return x*𝕦.m
    end
end

physical_units(x::Vector{L}, t::T) where {L<:Real, T<:Real} =
    physical_units.(x,:l), physical_units(t,:t)
physical_units(x::Vector{L}, v::Vector{V}) where {L<:Real, V<:Real} =
    physical_units.(x,:l), physical_units_v.(v,:v)
physical_units(x::Vector{L}, v::Vector{V}, t::T) where {L<:Real, V<:Real,T<:Real} =
    physical_units.(x,:l), physical_units.(v,:v), physical_units(t,:t)


"""Cosmos configuration"""
@with_kw struct CosmosConfig{F<:Unitful.Frequency, D<:Unitful.Density} <:AbstractConfig
    H₀::F = 67.66u"km/Mpc/s"
    ρ_c::D = uconvert(u"Msun/kpc^3", 3H₀^2/(8π*u"G") )
end


