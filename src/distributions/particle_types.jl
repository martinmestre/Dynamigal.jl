"""Particle types"""

@with_kw struct TestParticle <: AbstractTestParticle
    event::Event
end
TestParticle(t::T, x::Vector{D}, v::Vector{F}) where {T<:Real,D<:Real,F<:Real} =
    TestParticle(Event(t,x,v))
TestParticle(x::Vector{D}, v::Vector{F}) where {D<:Real,F<:Real} = TestParticle(Event(x,v))
TestParticle(t::T, x::Vector{D}, v::Vector{F}) where {T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} =
    TestParticle(Event(t,x,v))
TestParticle(x::Vector{D}, v::Vector{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    TestParticle(Event(x,v))


@with_kw struct Particle{T<:Real} <: AbstractParticle
    m::T = 1.0
    event::Event = Event()
end
Particle(m::M, t::T, x::Vector{D}, v::Vector{F}) where {M<:Real,T<:Real,D<:Real,F<:Real} =
    Particle(m, Event(t, x, v))
Particle(m::M, x::Vector{D}, v::Vector{F}) where {M<:Real, D<:Real,F<:Real} =
    Particle(m, Event(x, v))
Particle(x::Vector{D}, v::Vector{F}) where {D<:Real,F<:Real} =
    Particle(1.0, Event(x, v))
Particle(m::M, t::T, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} =
    Particle(ustrip(uconvert(𝕦.m, m)), Event(t, x, v))
Particle(m::M, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Velocity} =
    Particle(ustrip(uconvert(𝕦.m, m)), Event(x, v))


@with_kw struct MacroParticle{P<:AbstractPotential} <: AbstractMacroParticle
    pot::P
    event::Event = Event()
end
MacroParticle(p::P) where {P<:AbstractPotential} = MacroParticle(pot=p)
MacroParticle(p::P, t::T, x::Vector{D}, v::Vector{F}) where {P<:AbstractPotential,T<:Real,D<:Real,F<:Real} =
    MacroParticle(p, Event(t, x, v))
MacroParticle(p::P, x::Vector{D}, v::Vector{F}) where {P<:AbstractPotential, D<:Real,F<:Real} =
    MacroParticle(p, Event(x, v))
Particle(p::P, t::T, x::Vector{D}, v::Vector{F}) where {P<:AbstractPotential, T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} =
    MacroParticle(p, Event(t, x, v))
MacroParticle(p::P, x::Vector{D}, v::Vector{F}) where {P<:AbstractPotential, D<:Unitful.Length, F<:Unitful.Velocity} = MacroParticle(p, Event(x, v))



# """MacroParticleSystem"""
# mutable struct MacroParticleSystem{T, A} <: AbstractMacroParticle
#     macroparticles::T  # NTuple{N, AbstractMacroParticle}
#     accelerations::A

#     function MacroParticleSystem(macroparticles::Tuple{Vararg{AbstractMacroParticle}})
#         N = length(macroparticles)
#         if N == 1
#             throw(ArgumentError("MacroParticleSystem requires at least 2 elements, got $N"))
#         end

#         R = eltype(first(macroparticles).event.x)
#         accelerations = zeros(MVector{3*N, R})

#         return new{typeof(macroparticles), typeof(accelerations)}(macroparticles, accelerations)
#     end
# end

# MacroParticleSystem(macroparticles::AbstractMacroParticle...) =
#     MacroParticleSystem(macroparticles)

"""MacroParticleSystem"""
@with_kw mutable struct MacroParticleSystem{P <: NTuple{N, AbstractMacroParticle} where N, A} <: AbstractMacroParticle
    macroparticles::P
    accelerations::A
end

function MacroParticleSystem(macroparticles::NTuple{N, AbstractMacroParticle}) where {N}
    if N == 1
        throw(ArgumentError("MacroParticleSystem requires at least 2 elements, got $N"))
    end
    R = eltype(first(macroparticles).event.x)
    accelerations = zeros(MVector{3N, R})
    return MacroParticleSystem{typeof(macroparticles), typeof(accelerations)}(macroparticles, accelerations)
end

MacroParticleSystem(macroparticles::AbstractMacroParticle...) = MacroParticleSystem(macroparticles)
