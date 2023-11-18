module GalacticDynamics

using Reexport
using DifferentialEquations
using Parameters
@reexport using Zygote
using StaticArrays
@reexport using Unitful, UnitfulAstro


export SolverConfig
export potential, acceleration, ode, evolve
export Particle, TestParticle, Event
export Plummer, AllenSantillanHalo, MiyamotoNagaiDisk, PointMass
export example_Plummer, example_MiyamotoNagai, example_sum_of_potentials
export example_AllenSantillan


include("constants.jl")
include("config.jl")
include("abstract_types.jl")
include("potential_types.jl")
include("particle_types.jl")
include("orbit_types.jl")
include("potentials.jl")
include("accelerations.jl")
include("evolutions.jl")
include("examples.jl")
# include("metrics.jl")
# iclude("transformations.jl")
# include("reference_frames.jl")

end
