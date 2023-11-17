"""Orbit types"""


@with_kw struct Event{T<:U.Time, L<:U.Length, V<:U.Velocity} <: AbstractEvent
    t::T
    x::Vector{L}
    v::Vector{V}
end
Event(t::T, x::Vector{T}, v::Vector{T}) where {T<:Real} = Event(t*u_T, x*u_L, v*u_V)

@with_kw struct Orbit{T<:U.Time, L<:U.Length, V<:U.Velocity} <: AbstractOrbit
    t::Vector{T}
    x::Matrix{L}
    v::Matrix{V}
end
Orbit(t::vector{T}, x::Matrix{T}, v::Matrix{T}) where {T<:Real} = Orbit(t*u_T, x*u_L, v*u_V)

getindex(orb::Orbit, i) = Event(orb.t[i], orb.x[1:3,i], orb.v[1:3,i])
firstindex(orb::Orbit) = Event(orb.t[begin], orb.x[1:3,begin], orb.v[1:3,begin])
lastindex(orb::Orbit) = Event(orb.t[end], orb.x[1:3,end], orb.v[1:3,end])