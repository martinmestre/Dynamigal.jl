"""Types"""

import Base.:+

"""Abstract types"""
abstract type AbstractCosmos end

abstract type AbstractGeometry <: AbstractCosmos end
abstract type AbstractSpaceTime <: AbstractGeometry end
abstract type AbstractPotential <: AbstractGeometry end

abstract type AbstractDistribution <: AbstractCosmos end
abstract type AbstractDiscrete <: AbstractDistribution end
abstract type AbstractContinuous <: AbstractDistribution end
abstract type AbstractParticle <: AbstractDiscrete end
abstract type AbstractEnsemble <: AbstractDiscrete end
abstract type AbstractElementaryParticle <: AbstractParticle end
abstract type AbstractMacroParticle <: AbstractParticle end
abstract type AbstractStar <: AbstractMacroParticle end
abstract type AbstractDarkMatter <: AbstractMacroParticle end
abstract type AbstractPointGlobularCluster <: AbstractMacroParticle end
abstract type AbstractPointGalaxy <: AbstractMacroParticle end
abstract type AbstractGlobularCluster <: AbstractEnsemble end
abstract type AbstractGalaxy <: AbstractEnsemble end
abstract type AbstractContinuousGlobularCluster <: AbstractContinuous end
abstract type AbstractContinuousGalaxy <:AbstractContinuous end


Base.:+(a<:AbstractConservative, b<:AbstractConservative)::Vector{<:AbstractConservative} = [a,b]


