"""Evolution functions"""


"""Evolution of a an initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::D, v::F,
    Δt::Tuple{T,T}; options=SolverConfig()) where {D<:Real, F<:Real, T<:Real}
    (; solver, abstol, reltol ) = options
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, Δt, p)
    sol=solve(prob, solver; abstol=abstol, reltol=reltol)
    orb = Orbit(sol.t, sol[1:3,:], sol[4:6,:])
    return orb
end

"""Evolution of a unitful initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, v::Vector{<:Unitful.Velocity},
    Δt::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...)
    x = ustrip(uconvert.(𝕦.l, x))
    v = ustrip(uconvert.(𝕦.v, v))
    Δt = uconvert.(𝕦.τ, Δt)./𝕦.t
    return evolve(pot, x, v, Δt; kwargs...)
end

"""Evolution of an Event in an AbstractPotential"""
function evolve(pot::P, event::Event, Δt::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...) where {P<:UnionAbstractPotentials}
    Δt = uconvert.(𝕦.τ, Δt)./𝕦.t .+ event.t
    x = event.x
    v = event.v
    return evolve(pot, x, v, Δt, kwargs...)
end


"""Evolution of a TestParticle in an AbstractPotential"""
function evolve(pot::P, p::TestParticle, Δt::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...) where {P<:UnionAbstractPotentials}
    Δt = uconvert.(𝕦.τ, Δt)./𝕦.t .+ p.event.t
    x = p.event.x
    v = p.event.v
    return evolve(pot, x, v, Δt; kwargs...)
end

