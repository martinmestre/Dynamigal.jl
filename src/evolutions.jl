"""Evolution functions"""


"""Evolution of a an initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::Vector{D}, v::Vector{F},
   t_span::Tuple{T,T}; options=SolverConfig()) where {D<:Real, F<:Real, T<:Real}
    (; solver, abstol, reltol ) = options
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
    sol = solve(prob, solver; abstol=abstol, reltol=reltol)
    orb = Orbit(sol.t, sol[1:3,:], sol[4:6,:])
    return orb
end

"""Evolution of a unitful initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, v::Vector{<:Unitful.Velocity},
    t_span::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...)
    x, v = code_units(x, v)
    t_span = code_units.(t_span)
    return evolve(pot, x, v, t_span; kwargs...)
end

"""Evolution of an Event in an AbstractPotential"""
function evolve(pot::P, event::Event, t_span::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...) where {P<:UnionAbstractPotentials}
    t_span = code_units.(t_span) .+ event.t
    x = event.x
    v = event.v
    return evolve(pot, x, v, t_span; kwargs...)
end


"""Evolution of a TestParticle in an AbstractPotential"""
function evolve(pot::P, p::TestParticle, t_span::Tuple{<:Unitful.Time, <:Unitful.Time}; kwargs...) where {P<:UnionAbstractPotentials}
    t_span = code_units.(t_span) .+ p.event.t
    x = p.event.x
    v = p.event.v
    return evolve(pot, x, v, t_span; kwargs...)
end


"""Evolution of a system of MacroParticle"""
function evolve(mps::Vector{P}, t_span::Tuple{T,T}; options=SolverConfig()) where {P<:AbstractMacroParticle,T<:Real}
    (; solver, abstol, reltol, saveat) = options
    p = mps
    x = vcat([[mps[i].event.x for i ∈ eachindex(mps)]...]...)
    v = vcat([[mps[i].event.v for i ∈ eachindex(mps)]...]...)
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
    sol  =solve(prob, solver; abstol=abstol, reltol=reltol, saveat=saveat)
    sys_orb = Vector{Orbit}(undef, length(p))
    n = length(x)
    for i ∈ eachindex(p)
        j_x = selec(i)
        j_v = n+j_x
        sys_orb[i] = Orbit(sol.t, sol[j_x:j_x+2,:], sol[j_v:j_v+2,:])
    end
    return sys_orb
end