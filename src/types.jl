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
@with_kw mutable struct Plummer{M,L,V} <: AbstractPotential
        m::M
        b::L
        x::Vector{L}
        v::Vector{V}
end
function potential(pot::Plummer, x::AbstractArray)
    return -u.G * pot.m / sqrt(pot.b^2 + sum((x-pot.x).^2))
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

function acceleration(pot::AbstractPotential, x::AbstractArray)
    return -gradient(y->potential(pot, y), x)
end

function ode(u,p,t)
    pot = p[1]
    return SA[u[4:6]...,acceleration(pot, x[1:3])...]
end

function test()
    u₀ = SA[[1.0, 0.0, 0.0]ua.kpc...,[0.0,10.0,0.0]u"km/s"...]
    tspan = (0.0, 100.0)*u.Gyr
    prob = ODEProblem(ode, u₀, tspan)
    @btime sol=solve(prob, Tsit5())
    return sol
end