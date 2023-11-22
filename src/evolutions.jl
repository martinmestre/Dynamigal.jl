"""Evolution functions"""


"""Evolution of a an initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::D, v::F,
   t_span::Tuple{T,T}; options=SolverConfig()) where {D<:Real, F<:Real, T<:Real}
    (; solver, abstol, reltol ) = options
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀,t_span, p)
    sol=solve(prob, solver; abstol=abstol, reltol=reltol)
    orb = Orbit(sol.t, sol[1:3,:], sol[4:6,:])
    return orb
end

"""Evolution of a unitful initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, v::Vector{<:Unitful.Velocity},
   t_span::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...)
    x, v = code_units(x, v)
    t_span = code_units.(tspan)
    return evolve(pot, x, v, t_span; kwargs...)
end

"""Evolution of an Event in an AbstractPotential"""
function evolve(pot::P, event::Event, t_span::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...) where {P<:UnionAbstractPotentials}
    t_span = code_units.(t_span) .+ event.t
    x = event.x
    v = event.v
    return evolve(pot, x, v, t_span, kwargs...)
end


"""Evolution of a TestParticle in an AbstractPotential"""
function evolve(pot::P, p::TestParticle, t_span::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...) where {P<:UnionAbstractPotentials}
    t_span = code_units.(t_span) .+ p.event.t
    x = p.event.x
    v = p.event.v
    return evolve(pot, x, v, t_span; kwargs...)
end

