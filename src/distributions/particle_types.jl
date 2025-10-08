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




"""MacroParticleSystem"""

@with_kw mutable struct MacroParticleSystem{N, A} <: AbstractMacroParticleSystem
    macroparticles::NTuple{N,MacroParticle}
    accelerations::A
end

function MacroParticleSystem(macroparticles::NTuple{N, MacroParticle}) where {N}
    if N == 1
        throw(ArgumentError("MacroParticleSystem requires at least 2 elements, got $N"))
    end
    R = eltype(first(macroparticles).event.x)
    accelerations = zeros(MVector{3N, R})
    return MacroParticleSystem{N, typeof(accelerations)}(macroparticles, accelerations)
end
MacroParticleSystem(macroparticles::MacroParticle...) = MacroParticleSystem(macroparticles)

function MacroParticleSystem(galactic::P) where {P<:AbstractGalacticSystem}
    macroparticles = ntuple(i -> getfield(galactic, i), fieldcount(P))
    return MacroParticleSystem(macroparticles)
end


@with_kw struct LargeCloudMW{P,T} <: AbstractGalacticSystem
    mw::P
    cloud::T
end
function LargeCloudMW(mps::T) where {T<:AbstractMacroParticleSystem}
    return LargeCloudMW{typeof(mps[1]),typeof(mps[2]) }(mw=mps[1], cloud=mps[2])
end

@with_kw struct CloudsMW{P,T,W} <: AbstractGalacticSystem
    mw::P
    large::T
    small::W
end
function CloudsMW(mps::T) where {T<:AbstractMacroParticleSystem}
    return CloudsMW{typeof(mps[1]),typeof(mps[2]),typeof(mps[3])}(mw=mps[1], large=mps[2], small=mps[3])
end

@with_kw struct SagCloudsMW{P,T,W,Q} <: AbstractGalacticSystem
    mw::P
    large::T
    small::W
    sag::Q
end
function SagCloudsMW(mps::T) where {T<:AbstractMacroParticleSystem}
    return SagCloudsMW{typeof(mps[1]),typeof(mps[2]),typeof(mps[3]),typeof(mps[4])}(mw=mps[1], large=mps[2], small=mps[3], sag=mps[4])
end