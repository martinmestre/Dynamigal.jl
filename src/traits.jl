
"""Traits"""


abstract type SystemTrait end
struct GenSys <: SystemTrait end
struct GenPerfSys <: SystemTrait end
struct CloudsMW <: SystemTrait end
struct SagCloudsMW <: SystemTrait end

SystemTrait(::Type) = GenSys()
# then in code I do if and only if I want to use the customize function:
# SystemTrait(::Type{typeof(mps)}) = CloudsMW()
# warning: this sets the same SystemTrait function for all the struct of the same
# type of "mps".

"""Macro to set system trait"""
macro set_trait(system, trait)
    return :(SystemTrait(::Type{$(typeof(system))}) = $trait())
end

# luego se corre:
# Más simple y readable
# @set_trait system CloudsMW

# Para cambiar
# @set_trait system SagCloudsMW


