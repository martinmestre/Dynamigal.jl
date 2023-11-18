"""Orbit types"""

@with_kw struct Event{D<:Real, T::typeof(u_T), L::typeof(u_L), V::typeof(u_V)} <: AbstractEvent
    t::D
    x::Vector{D}
    v::Vector{D}
    units::Tuple{T,L,V} = (u_T, u_L, u_V)
end


# @with_kw struct Event{T<:U.Time, L<:U.Length, V<:U.Velocity, D<:Real} <: AbstractEvent
#     t_u::T
#     x_u::Vector{L}
#     v_u::Vector{V}
#     t::D
#     x::Vector{D}
#     v::Vector{D}
#     function Event{T,L,V}(t_u, x_u, v_u) where {T,L,V}
#         t = t_u / uconvert(unit(t_u), u_T)
#         x = ustrip(uconvert.(u_L, x_u))
#         v = ustrip(uconvert.(u_V, v_u))
#         D = typeof(t)
#         return new{T,L,V,D}(t_u, x_u, v_u, t, x, v)
#     end
# end
# Event(t_u::T, x_u::Vector{L}, v_u::Vector{V}) where {T,L,V} = Event{T,L,V}(t_u, x_u, v_u)
# Event(t::T, x::Vector{T}, v::Vector{T}) where {T<:Real} = Event(t*u_T, x*u_L, v*u_V)


@with_kw struct Orbit{T<:U.Time, L<:U.Length, V<:U.Velocity} <: AbstractOrbit
    t::Vector{T}
    x::Matrix{L}
    v::Matrix{V}
end
Orbit(t::Vector{T}, x::Matrix{T}, v::Matrix{T}) where {T<:Real} = Orbit(t*u_T, x*u_L, v*u_V)

getindex(orb::Orbit, i) = Event(orb.t[i], orb.x[1:3,i], orb.v[1:3,i])
firstindex(orb::Orbit) = Event(orb.t[begin], orb.x[1:3,begin], orb.v[1:3,begin])
lastindex(orb::Orbit) = Event(orb.t[end], orb.x[1:3,end], orb.v[1:3,end])