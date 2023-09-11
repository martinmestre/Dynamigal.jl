"""Types"""

import Base.:+

"""Abstract types"""
abstract type AbstractCosmos end

abstract type AbstractGeometry <: AbstractCosmos end
abstract type AbstractSpaceTime <: AbstractGeometry end
abstract type AbstractPotential <: AbstractGeometry end

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


# Base.:+(a::Union{<:AbstractPotential,Vector{<:AbstractPotentail}},
#         b::Union{<:AbstractPotential,Vector{<:AbstractPotentail}})
#         ::Vector{<:AbstractPotential} = vcat(a,b)


"""Concrete types (structs)"""

@with_kw mutable struct Plummer{T<:Real} <: AbstractPotential
        m_u::Quantity{T}
        b_u::Quantity{T}
        m::T
        b::T
        function Plummer{T}(m_u, b_u) where {T}
            m = uconvert(u"Msun", m_u).val
            b = uconvert(u"kpc", b_u).val
            return new{T}(m_u, b_u, m, b)
        end
end
Plummer(m_u::Quantity{T}, b_u::Quantity{T}) where {T} = Plummer{T}(m_u, b_u)

function potential(pot::Plummer, x::Vector{<:Quantity{T}}) where {T}
    return -u"G"*pot.m_u / sqrt(pot.b_u^2 + x'x)
end

function potential(pot::Plummer, x::AbstractArray{T}) where {T}
    return -G*pot.m / sqrt(pot.b^2 + x'x)
end

@with_kw mutable struct Particle{M,L,V} <: AbstractMacroParticle
        m::M
        x::Vector{L}
        v::Vector{V}
end

@with_kw mutable struct TestParticle{L,V} <: AbstractTestParticle
        x::Vector{L}
        v::Vector{V}
end

function acceleration(pot::AbstractPotential, x::Vector{<:Quantity{T}}) where {T}
    x = [uconvert(u"kpc", x[i]).val for i ∈ eachindex(x)]
    return -1.0.*gradient(y->potential(pot, y), x)[1]
end
function acceleration(pot::AbstractPotential, x::AbstractArray{T}) where {T}
    return -1.0.*gradient(y->potential(pot, y), x)[1]
end

function ode(u,p,t)
    pot = p[1]
    return SA[u[4:6]...,acceleration(pot, u[1:3])...]
end

function test()
    plum = Plummer(10.0^8*u"Msun",10.0u"kpc")
    p = [plum]
    x₀ = [10.0, 0.0, 0.0]
    v₀ = [0.0,10.0,0.0]
    u₀ = SA[x₀...,v₀...]
    tspan = (0.0, 100.0)
    prob = ODEProblem(ode, u₀, tspan, p)
    sol=solve(prob, Tsit5(); abstol=5.e-8, reltol=5.e-8)
    return sol
end