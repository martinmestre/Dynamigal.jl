module GalacticDynamics

using Reexport
using DifferentialEquations
using Parameters
@reexport using Zygote
using StaticArrays
@reexport using Unitful, UnitfulAstro
import UnitfulChainRules
U=Unitful

export G, u_T
export Plummer, potential, acceleration
export evolve
export Particle, TestParticle
export ode, test

include("constants.jl")
include("types.jl")
include("methods.jl")
# include("metrics.jl")
# include("potentials.jl")
# iclude("transformations.jl")
# include("reference_frames")

end
