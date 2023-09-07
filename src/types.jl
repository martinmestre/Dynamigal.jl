"""Types"""

import Base.:+

"""Abstract types"""
abstract type AbstractCosmos end
abstract type AbstractSpaceTime <: AbstractCosmos end
abstract type AbstractConservative<: AbstractSpaceTime end
abstract type AbstractNonConservative <: AbstractSpaceTime end
abstract type AbstractPotential <: AbstractConservative end

Base.:+(a<:AbstractConservative, b<:AbstractConservative)::Vector{<:AbstractConservative} = [a,b]


