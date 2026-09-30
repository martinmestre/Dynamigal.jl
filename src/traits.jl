
"""Traits"""

"""Evolution/ODE traits"""
abstract type SystemTrait end
struct GenSysTrait <: SystemTrait end
struct GenSysMutOdeTrait <: SystemTrait end
struct RawSolutionTrait <: SystemTrait end
struct GalacticTrait <: SystemTrait end
struct PerfGalacticTrait <: SystemTrait end
struct MutualFrictionTrait <: SystemTrait end

SystemTrait(::Type) = GenSysTrait()
SystemTrait(::Type{<:AbstractGalacticSystem}) = GalacticTrait()

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
# @set_trait system GenSysMutOdeTrait

# Para cambiar
# @set_trait system GenSysTrait

"""For specifying a not-default derivation/integration of a quantity"""
abstract type FromTrait end
struct FromPotentialTrait <: FromTrait end
struct FromMassTrait <: FromTrait end
struct FromAccelerationTrait <: FromTrait end
struct FromDensityTrait <: FromTrait end
struct FromDistributionTrait <: FromTrait end