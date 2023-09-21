module GalacticDynamics

using Reexport
using DifferentialEquations
using Parameters
@reexport using Zygote
using StaticArrays
@reexport using Unitful, UnitfulAstro
import UnitfulChainRules
const U=Unitful

export G, u_T
export Plummer, potential, acceleration, ode, evolve
export Particle, TestParticle
export example_Plummer, example_MiyamotoNagai, example_sum_of_potentials
export example_AllenSantillan

include("constants.jl")
include("types.jl")
include("potentials.jl")
include("acceleration.jl")
include("evolutions.jl")
include("examples.jl")
# include("metrics.jl")
# include("potentials.jl")
# iclude("transformations.jl")
# include("reference_frames")

end
