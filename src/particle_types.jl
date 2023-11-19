"""Particle types"""


@with_kw struct TestParticle <: AbstractTestParticle
    event::Event
end
TestParticle(t::T, x::Vector{D}, v::Vector{F}) where {T<:Real,D<:Real,F<:Real} =
    TestParticle(Event(t,x,v))
TestParticle(x::Vector{D}, v::Vector{F}) where {D<:Real,F<:Real} = TestParticle(Event(x,v))
TestParticle(t::T, x::Vector{D}, v::Vector{F}) where {T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} = TestParticle(Event(t,x,v))
TestParticle(x::Vector{D}, v::Vector{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    TestParticle(Event(x,v))


@with_kw struct Particle{T<:Real} <: AbstractParticle
    m::T = 1.0
    event::Event
end
Particle(m::M, t::T, x::Vector{D}, v::Vector{F}) where {M<:Real,T<:Real,D<:Real,F<:Real} =
    Particle(m, Event(t, x, v))
Particle(m::M, t::T, x::Vector{D}, v::Vector{F}) where {M<:Real,T<:Real,D<:Real,F<:Real} =
    Particle(m, Event(t, x, v))
Particle(m::M, x::Vector{D}, v::Vector{F}) where {M<:Real, D<:Real,F<:Real} =
    Particle(m, Event(x, v))
Particle(m::M, t::T, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} = Particle(m, Event(t, x, v))
Particle(m::M, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Velocity} =
    Particle(m, Event(x, v))



