"""Types"""

import Base.:+

"""Abstract types"""
abstract type AbstractCosmos end

abstract type AbstractGeometry <: AbstractCosmos end
abstract type AbstractSpaceTime <: AbstractGeometry end
abstract type AbstractPotential <: AbstractGeometry end
abstract type AbstractDiskPotential <: AbstractPotential end
abstract type AbstractBulgePotential <: AbstractPotential end
abstract type AbstractHaloPotential <: AbstractPotential end


abstract type AbstractDistribution <: AbstractCosmos end
abstract type AbstractDiscreteDistribution <: AbstractDistribution end
abstract type AbstractContinuousDistribution <: AbstractDistribution end

abstract type AbstractParticle <: AbstractDiscreteDistribution end
abstract type AbstractEnsemble <: AbstractDiscreteDistribution end
abstract type AbstractElementaryParticle <: AbstractParticle end
abstract type AbstractTestParticle <: AbstractParticle end
abstract type AbstractMacroParticle <: AbstractParticle end

abstract type AbstractGlobularCluster <: AbstractEnsemble end
abstract type AbstractGalaxy <: AbstractEnsemble end
abstract type AbstractHalo <: AbstractEnsemble end
abstract type AbstractDisk <: AbstractEnsemble end
abstract type AbstractBulge <: AbstractEnsemble end

abstract type AbstractContinuousGlobularCluster <: AbstractContinuousDistribution end
abstract type AbstractContinuousGalaxy <:AbstractContinuousDistribution end
abstract type AbstractContinuousHalo <:AbstractContinuousDistribution end
abstract type AbstractContinuousSubHalo <:AbstractContinuousDistribution end
abstract type AbstractContinuousDisk <:AbstractContinuousDistribution end
abstract type AbstractContinuousBulge <:AbstractContinuousDistribution end

abstract type AbstractOrbit <: AbstractCosmos end

const UnionAbstractPotentials = Union{AbstractPotential, Vector{<:AbstractPotential}}

"""Overloading sum in order to sum potentials"""
function Base.:+(a::UnionAbstractPotentials, b::UnionAbstractPotentials)
    return vcat(a,b)
end

"""Concrete types (structs)"""

@with_kw struct Plummer{M<:U.Mass, L<:U.Length,T<:Real} <: AbstractPotential
        m_u::M
        b_u::L
        m::T
        b::T
        function Plummer{M,L}(m_u, b_u) where {M,L}
            m = uconvert(u_M, m_u).val
            b = uconvert(u_L, b_u).val
            T=typeof(m)
            return new{M,L,T}(m_u, b_u, m, b)
        end
end
Plummer(m_u::M, b_u::L) where {M,L} = Plummer{M,L}(m_u, b_u)

@with_kw struct MiyamotoNagaiDisk{M<:U.Mass, L<:U.Length,T<:Real} <: AbstractDiskPotential
    m_u::M
    a_u::L
    b_u::L
    m::T
    a::T
    b::T
    function MiyamotoNagaiDisk{M,L}(m_u, a_u, b_u) where {M,L}
        m = uconvert(u_M, m_u).val
        a = uconvert(u_L, a_u).val
        b = uconvert(u_L, b_u).val
        T=typeof(m)
        return new{M,L,T}(m_u, a_u, b_u, m, a, b)
    end
end
MiyamotoNagaiDisk(m_u::M, a_u::L, b_u::L) where {M,L} = MiyamotoNagaiDisk{M,L}(m_u, a_u, b_u)


@with_kw struct AllenSantillanHalo{M<:U.Mass, L<:U.Length, D<:Real, T<:Real} <: AbstractDiskPotential
    m_u::M
    a_u::L
    Λ_u::L
    γ::D
    m::T
    a::T
    Λ::T
    function AllenSantillanHalo{M,L,D}(m_u, a_u, Λ_u, γ) where {M,L,D}
        m = uconvert(u_M, m_u).val
        a = uconvert(u_L, a_u).val
        Λ = uconvert(u_L, Λ_u).val
        T=typeof(m)
        return new{M,L,D,T}(m_u, a_u, Λ_u, γ, m, a, Λ)
    end
end
AllenSantillanHalo(m_u::M, a_u::L, Λ_u::L, γ::D) where {M,L,D} = AllenSantillanHalo{M,L,D}(m_u, a_u, Λ_u, γ)


@with_kw struct Particle{M<:U.Mass,L<:U.Length,V<:U.Velocity,T<:Real} <: AbstractMacroParticle
        m_u::M
        x_u::Vector{L}
        v_u::Vector{V}
        m::T
        x::Vector{T}
        v::Vector{T}
        function Particle{M,L,V}(m_u, x_u, v_u) where {M,L,V}
            m = uconvert(u_M, m_u).val
            x = ustrip(uconvert.(u_L, x_u))
            v = ustrip(uconvert.(u_V, v_u))
            T=typeof(m)
            return new{M,L,V,T}(m_u, x_u, v_u, m, x, v)
        end
end
Particle(m_u::M, x_u::Vector{L}, v_u::Vector{V}) where {M,L,V} = Particle{M,L,V}(m_u, x_u, v_u)

@with_kw struct TestParticle{L<:U.Length,V<:U.Velocity,T<:Real} <: AbstractTestParticle
    x_u::Vector{L}
    v_u::Vector{V}
    x::Vector{T}
    v::Vector{T}
    function TestParticle{L,V}(x_u, v_u) where {L,V}
        x = ustrip(uconvert.(u_L, x_u))
        v = ustrip(uconvert.(u_V, v_u))
        T=typeof(x[begin])
        return new{L,V,T}(x_u, v_u, x, v)
    end
end
TestParticle(x_u::Vector{L}, v_u::Vector{V}) where {L,V} = TestParticle{L,V}(x_u, v_u)

@with_kw struct Orbit{T<:U.Time, L<:U.Length, V<:U.Velocity} <: AbstractOrbit
    t::Vector{T}
    x::Matrix{L}
    v::Matrix{V}
end
