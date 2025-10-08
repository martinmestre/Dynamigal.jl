
"""Traits"""

"""Evolution/ODE traits"""
abstract type SystemTrait end
struct GenSysTrait <: SystemTrait end
struct GenSysMutOdeTrait <: SystemTrait end

SystemTrait(::Type) = GenSysTrait()
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

"""Acceleration traits"""
abstract type FrictionTrait end
struct FrictionIncludedTrait <: FrictionTrait end
struct FrictionlessTrait <: FrictionTrait end

FrictionTrait(::Type) = FrictionlessTrait()

macro set_friction_trait(system, trait)
    return :(FrictionTrait(::Type{typeof($(esc(system)))}) = $(esc(trait))())
end
