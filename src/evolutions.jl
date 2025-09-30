"""Evolution functions"""


"""Evolution of a an initial condition in an AbstractPotential"""
function evolve(pot::P, x::AbstractVector{D}, v::AbstractVector{F},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential, D, F, T}
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
    sol = solve(prob, solver; options...)
    orb = Orbit(sol.t, sol[sis,:], sol[siss,:])
    return orb
end

function _evolve(pot::P, x::AbstractVector{D}, v::AbstractVector{F},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential, D, F, T}
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(_ode, u₀, t_span, p)
    sol = solve(prob, solver; options...)
    orb = Orbit(sol.t, sol[sis,:], sol[siss,:])
    return orb
end

"""Evolution of a unitful initial condition in an AbstractPotential"""
function evolve(pot::P, x::Vector{<:Unitful.Length}, v::Vector{<:Unitful.Velocity},
    t_span::Tuple{<:Unitful.Time, <:Unitful.Time}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential}
    x, v = adimensional(x, v)
    t_span = adimensional.(t_span)
    return evolve(pot, x, v, t_span, solver; options)
end

"""Evolution of an Event in an AbstractPotential"""
function evolve(pot::P, event::Event, t_span::Tuple{<:Unitful.Time, <:Unitful.Time}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential}
    t_span = adimensional.(t_span) .+ event.t
    x = event.x
    v = event.v
    return evolve(pot, x, v, t_span, solver; options)
end


"""Evolution of a TestParticle in an AbstractPotential"""
function evolve(pot::P, p::TestParticle, t_span::Tuple{<:Unitful.Time, <:Unitful.Time}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential}
    t_span = code_units.(t_span) .+ p.event.t
    x = p.event.x
    v = p.event.v
    return evolve(pot, x, v, t_span, solver; options)
end


"""Evolution of a MacroParticleSystem, main method"""
function evolve(mps::T, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {T,R<:Real}
    @show T SystemTrait(T)
    return evolve(SystemTrait(T), mps, t_span, solver; options=options)
end

"""Evolution of a system of MacroParticle, general type"""
function evolve(::GenSys, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    println("Caso GenSys")
    @show options
    x = vcat([[mps[i].event.x for i ∈ eachindex(mps)]...]...)
    v = vcat([[mps[i].event.v for i ∈ eachindex(mps)]...]...)
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, mps)
    sol  = solve(prob, solver; options...)
    sys_orb = Vector{Orbit}(undef, length(mps))
    n = length(x)
    for i ∈ eachindex(mps)
        j_x = selec(i)
        j_v = n+j_x
        sys_orb[i] = Orbit(sol.t, sol[j_x:j_x+2,:], sol[j_v:j_v+2,:])
    end
    return sys_orb
end

"""Evolution of a system of MacroParticle, general performant type"""
function evolve(::GenPerfSys, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    println("Caso GenPerfSys")
    @show options
    x = vcat([[mps[i].event.x for i ∈ eachindex(mps)]...]...)
    v = vcat([[mps[i].event.v for i ∈ eachindex(mps)]...]...)
    u₀ = [x...,v...]
    prob = ODEProblem(ode!, u₀, t_span, mps)
    sol  =solve(prob, solver; options...)
    sys_orb = Vector{Orbit}(undef, length(mps))
    n = length(x)
    for i ∈ eachindex(mps)
        j_x = selec(i)
        j_v = n+j_x
        sys_orb[i] = Orbit(sol.t, sol[j_x:j_x+2,:], sol[j_v:j_v+2,:])
    end
    return sys_orb
end

"""Evolution of a system of MacroParticle, Clouds & MW type"""
#evolve(::CloudsMW, mps::T, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {T,R<:Real} =




