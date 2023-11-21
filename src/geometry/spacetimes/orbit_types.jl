"""Orbit types"""


@with_kw struct Event{T<:Real,D<:Real,F<:Real} <: AbstractEvent
    t::T = 0.0
    x::Vector{D} = zeros(3)
    v::Vector{F} = zeros(3)
end
Event(x::Vector{D}, v::Vector{F}) where {D<:Real,F<:Real} = Event(0.0, x, v)
Event(t::T, x::Vector{D}, v::Vector{F}) where {T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} =
    Event( ustrip(uconvert(𝕦.t, t)),  ustrip(uconvert.(𝕦.l, x)), ustrip(uconvert.(𝕦.v, v)) )
Event(x::Vector{D}, v::Vector{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    Event( ustrip(uconvert.(𝕦.l, x)), ustrip(uconvert.(𝕦.v, v)) )


    @with_kw struct Orbit{T<:Real,D<:Real,F<:Real} <: AbstractOrbit
    t::Vector{T}
    x::Matrix{D}
    v::Matrix{F}
end
Orbit(t::Vector{T}, x::Matrix{D}, v::Matrix{F}) where {T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} = Orbit( ustrip(uconvert(𝕦.t, t)),  ustrip(uconvert.(𝕦.l, x)), ustrip(uconvert.(𝕦.v, v)) )

getindex(orb::Orbit, i) = Event(orb.t[i], orb.x[1:3,i], orb.v[1:3,i])
firstindex(orb::Orbit) = Event(orb.t[begin], orb.x[1:3,begin], orb.v[1:3,begin])
lastindex(orb::Orbit) = Event(orb.t[end], orb.x[1:3,end], orb.v[1:3,end])