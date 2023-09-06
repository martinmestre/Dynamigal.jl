"""Types"""

import Base.:+

"""Abstract types"""
abstract type AbstractCosmos end
abstract type AbstractSpaceTime <: AbstractCosmos end
abstract type AbstractConservativeForce <: AbstractSpaceTime end
abstract type AbstractNonConservativeForce <: AbstractSpaceTime end
abstract type AbstractPotential <: AbstractConservativeForce end

Base.:+(a<:AbstractConservativeForce, b<:AbstractConservativeForce)::Vector{<:AbstractConservativeForce} = [a,b]


