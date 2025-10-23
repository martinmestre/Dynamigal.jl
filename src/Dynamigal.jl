module Dynamigal

using Reexport
@reexport using OrdinaryDiffEq
using Parameters
@reexport using Zygote
@reexport using StaticArrays
@reexport using Unitful, UnitfulAstro
using NamedTupleTools
@reexport using LinearAlgebra
@reexport using SpecialFunctions

export SolverConfig, UnitsConfig, CosmosConfig, SolverOptions, FrictionConfig
export ntSolverOptions
export 𝕦, G, 𝕤, 𝕔, 𝕗, H₀, sis, six
export potential
export acceleration, acceleration!, acceleration_c!
export ode, _ode, ode_c, ode_perf
export evolve, _evolve, evolve_c
export density
export circular_velocity
export Particle, TestParticle, MacroParticle
export MacroParticleSystem
export LargeCloudMW, CloudsMW, SagCloudsMW
export Event, Orbit, Snapshot
export TimeDependent, Kepler, Plummer, AllenSantillanHalo, MiyamotoNagaiDisk
export Hernquist, NFW
export OscillatoryKepler
export CompositePotential
export concentration
export example_Plummer, example_MiyamotoNagai, example_sum_of_potentials
export example_AllenSantillan
export example_of_mps
export example_cloudsMW
export example_cloudsMW_friction
export code_units, physical_units, adimensional
export r_vir_nfw
export ChandrasekharFriction
export drag
export MilkyWayBovy2014
export SystemTrait, @set_system_trait
export GenSysTrait, GenSysMutOdeTrait
export GalacticTrait, PerfGalacticTrait




include("abstract_types.jl")
include("traits.jl")
include("config.jl")
include("constants.jl")
include("geometry/potentials/potential_types.jl")
include("geometry/potentials/potentials.jl")
include("geometry/potentials/customized.jl")
include("geometry/spacetimes/orbit_types.jl")
include("distributions/particle_types.jl")
include("distributions/ensemble_types.jl")
include("dissipative/dissipative_types.jl")
include("dissipative/drags.jl")
include("odes.jl")
include("accelerations.jl")
include("evolutions.jl")
include("densities.jl")
include("circular_velocity.jl")
include("overloads.jl")
include("examples.jl")

# include("metrics.jl")
# iclude("transformations.jl")
# include("reference_frames.jl")

end
