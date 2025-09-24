"""Types"""


"""Abstract types"""
abstract type AbstractConfig end
abstract type AbstractCosmos end

abstract type AbstractGeometry <: AbstractCosmos end

abstract type AbstractPotential <: AbstractGeometry end
abstract type AbstractStaticPotential <: AbstractPotential end

abstract type AbstractSpaceTime <: AbstractGeometry end

abstract type AbstractOrbit <: AbstractSpaceTime end
abstract type AbstractEvent <: AbstractSpaceTime end

abstract type AbstractDistribution <: AbstractCosmos end

abstract type AbstractDiscreteDistribution <: AbstractDistribution end

abstract type AbstractParticle <: AbstractDiscreteDistribution end

abstract type AbstractElementaryParticle <: AbstractParticle end
abstract type AbstractTestParticle <: AbstractParticle end
abstract type AbstractMacroParticle <: AbstractParticle end

abstract type AbstractEnsemble <: AbstractDiscreteDistribution end

abstract type AbstractGlobularCluster <: AbstractEnsemble end
abstract type AbstractGalaxy <: AbstractEnsemble end


abstract type AbstractContinuousDistribution <: AbstractDistribution end

abstract type AbstractContinuousGlobularCluster <: AbstractContinuousDistribution end
abstract type AbstractContinuousGalaxy <:AbstractContinuousDistribution end

