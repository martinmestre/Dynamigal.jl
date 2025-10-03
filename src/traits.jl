
"""Traits"""


abstract type SystemTrait end
struct GenSys <: SystemTrait end
struct GenSysMutODE <: SystemTrait end

@with_kw struct LargeCloudMW{P<:AbstractMacroParticle} <: SystemTrait
    mw::P
    large::P
end
function LargeCloudMW(mps::T) where T<:AbstractMacroParticleSystem
    return LargeCloudMW{typeof(mps[1])}(mw=mps[1], large=mps[2])
end

@with_kw struct CloudsMW{P<:AbstractMacroParticle} <: SystemTrait
    mw::P
    large::P
    small::P
end
function CloudsMW(mps::T) where {T<:AbstractMacroParticleSystem}
    return CloudsMW{eltype(T)}(mw=mps[1], large=mps[2], small=mps[3])
end

@with_kw struct SagCloudsMW{P<:AbstractMacroParticle} <: SystemTrait
    mw::P
    large::P
    small::P
    sag::P
end
function SagCloudsMW(mps::T) where {T<:AbstractMacroParticleSystem}
    return SagCloudsMW{eltype(T)}(mw=mps[1], large=mps[2], small=mps[3], sag=mps[4])
end

SystemTrait(::Type) = GenSys()
# then in code I do if and only if I want to use the customize function:
# SystemTrait(::Type{typeof(mps)}) = CloudsMW()
# warning: this sets the same SystemTrait function for all the struct of the same
# type of "mps".

"""Macro to set system trait"""
macro set_system_trait(system, trait)
    return :(SystemTrait(::Type{typeof($(esc(system)))}) = $(esc(trait))())
end
# luego se corre:
# Más simple y readable
# @set_trait system CloudsMW

# Para cambiar
# @set_trait system SagCloudsMW


