module GalacticDynamics

using Reexport
using DifferentialEquations
using Parameters
@reexport using Zygote
using StaticArrays
@reexport using Unitful, UnitfulAstro
using NamedTupleTools
using BenchmarkTools


export SolverConfig, UnitsConfig, CosmosConfig, SolverOptions
export ntSolverOptions
export 𝕦, G, 𝕤, 𝕔, H₀
export potential, acceleration, ode, evolve
export circular_velocity
export Particle, TestParticle, MacroParticle
export Event, Orbit, Snapshot
export TimeDependent, Kepler, Plummer, AllenSantillanHalo, MiyamotoNagaiDisk
export Hernquist, NFW
export concentration
export example_Plummer, example_MiyamotoNagai, example_sum_of_potentials
export example_AllenSantillan
export code_units, physical_units, adimensional



include("config.jl")
include("constants.jl")
include("abstract_types.jl")
include("geometry/potentials/potential_types.jl")
include("geometry/potentials/potentials.jl")
include("geometry/spacetimes/orbit_types.jl")
include("distributions/particle_types.jl")
include("distributions/ensemble_types.jl")
include("odes.jl")
include("accelerations.jl")
include("evolutions.jl")
include("circular_velocity.jl")
include("examples.jl")

# include("metrics.jl")
# iclude("transformations.jl")
# include("reference_frames.jl")

end
