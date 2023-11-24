module GalacticDynamics

using Reexport
using DifferentialEquations
using Parameters
@reexport using Zygote
using StaticArrays
@reexport using Unitful, UnitfulAstro


export SolverConfig, UnitsConfig, 𝕦
export potential, acceleration, ode, evolve
export Particle, TestParticle
export Event, Orbit, Snapshot
export TimeDependent, PointMass, Plummer, AllenSantillanHalo, MiyamotoNagaiDisk
export example_Plummer, example_MiyamotoNagai, example_sum_of_potentials
export example_AllenSantillan



include("config.jl")
include("constants.jl")
include("abstract_types.jl")
include("geometry/potentials/potential_types.jl")
include("geometry/potentials/potentials.jl")
include("geometry/spacetimes/orbit_types.jl")
include("distributions/particle_types.jl")
include("distributions/ensemble_types.jl")
include("accelerations.jl")
include("evolutions.jl")
include("examples.jl")
# include("metrics.jl")
# iclude("transformations.jl")
# include("reference_frames.jl")

end
