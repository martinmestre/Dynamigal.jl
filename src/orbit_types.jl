"""Orbit types"""


@with_kw struct Orbit{T<:U.Time, L<:U.Length, V<:U.Velocity} <: AbstractOrbit
    t::Vector{T}
    x::Matrix{L}
    v::Matrix{V}
end

@with_kw struct PhaseSpacePoint{T<:U.Time, L<:U.Length, V<:U.Velocity} <: AbstractPhaseSpacePoint
    t::T
    x::Vector{L}
    v::Vector{V}
end
getindex(orb::Orbit,i) = PhaseSpacePoint(orb.t[i], orb.x[1:3,i], orb.v[1:3,i])
firstindex(orb::Orbit) = PhaseSpacePoint(orb.t[begin], orb.x[1:3,begin], orb.v[1:3,begin])
lastindex(orb::Orbit) = PhaseSpacePoint(orb.t[end], orb.x[1:3,end], orb.v[1:3,end])