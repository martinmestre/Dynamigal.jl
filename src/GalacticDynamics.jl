module GalacticDynamics

using DifferentialEquations
using Parameters
using Zygote
using StaticArrays
using Unitful, UnitfulAstro
import UnitfulChainRules
U=Unitful

export G
export Plummer, potential, acceleration
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
