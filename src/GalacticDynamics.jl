module GalacticDynamics

using DifferentialEquations
using Parameters
using Zygote
using StaticArrays
using Unitful, UnitfulAstro
import UnitfulChainRules
u=UnitfulAstro

export G
export Plummer, potential, acceleration
export ode, test
include("constants.jl")
include("types.jl")
# include("metrics.jl")
# include("potentials.jl")
# iclude("transformations.jl")
# include("reference_frames")

end
