module GalacticDynamics

using DifferentialEquations
using Parameters
using Zygote
using StaticArrays
using Unitful, UnitfulAstro
import UnitfulChainRules
u=Unitful
ua=UnitfulAstro

export Plummer, potential, acceleration
include("types.jl")
# include("metrics.jl")
# include("potentials.jl")
# iclude("transformations.jl")
# include("reference_frames")

end
