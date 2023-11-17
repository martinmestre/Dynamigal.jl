"""Types"""

import Base.:+
import Base.:getindex
import Base.:firstindex, Base.:lastindex

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
abstract type AbstractEvent <: AbstractCosmos end

const UnionAbstractPotentials = Union{AbstractPotential, Vector{<:AbstractPotential}}

"""Overloading sum in order to sum potentials"""
function Base.:+(a::UnionAbstractPotentials, b::UnionAbstractPotentials)
    return vcat(a,b)
end






