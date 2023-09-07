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
abstract type AbstractMacroParticle <: AbstractParticle end
abstract type AbstractStar <: AbstractMacroParticle end
abstract type AbstractDarkMatter <: AbstractMacroParticle end
abstract type AbstractPointGlobularCluster <: AbstractMacroParticle end
abstract type AbstractPointGalaxy <: AbstractMacroParticle end
abstract type AbstractPointSubHalo <: AbstractMacroParticle end
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


Base.:+(a<:AbstractPotential, b<:AbstractPotential)::Vector{<:AbstractPotential} = [a,b]
Base.:+(a::Vector{<:AbstractPotentail}, b<:AbstractPotential)::Vector{<:AbstractPotential} = vcat(a,b)
Base.:+(a::Vector{<:AbstractPotentail}, b::Vector{<:AbstractPotential})::Vector{<:AbstractPotential} = vcat(a,b)
