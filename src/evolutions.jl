"""Evolution functions"""


"""Evolution of a an initial condition in an AbstractPotential"""
function _evolve(pot::P, x::AbstractVector{D}, v::AbstractVector{F},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential, D, F, T}
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(_ode, u₀, t_span, p)
    sol = solve(prob, solver; options...)
    orb = Orbit(sol.t, sol[sis,:], sol[siss,:])
    return orb
end
# evolve(below) is a little faster than _evolve(above)
function evolve(pot::P, x::AbstractVector{D}, v::AbstractVector{F},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential, D, F, T}
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
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


"""Evolution of a P <: AbstractMacroParticleSystem, main method"""
function evolve(mps::T, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {T<:AbstractMacroParticleSystem, R<:Real}
    return evolve(SystemTrait(T), mps, t_span, solver; options=options)
end

"""Evolution of a system of MacroParticle, general type"""
function evolve(::GenSysTrait, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
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


"""Evolution of a system of MacroParticle using Mutating ODE trait"""
function evolve(::GenSysMutOdeTrait, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
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

"""Evolution of a LargeCloudMW (<: GalacticSystem)"""
function evolve(::GalacticTrait, cloudMW::LargeCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    x_mw = cloudMW.mw.event.x
    x_cl = cloudMW.cloud.event.x
    v_mw = cloudMW.mw.event.v
    v_cl = cloudMW.cloud.event.v
    u₀ = SVector{12,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3], x_cl[1], x_cl[2], x_cl[3],
                                    v_mw[1], v_mw[2], v_mw[3], v_cl[1], v_cl[2], v_cl[3])
    prob = ODEProblem(ode, u₀, t_span, cloudMW)
    sol  = solve(prob, solver; options...)
    sys_orb = Vector{Orbit}(undef, 2)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[7:9,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[10:12,:])
    return sys_orb
end

"""Evolution of a LargeCloudMW (<: GalacticSystem)"""
function evolve(::PerfGalacticTrait, cloudMW::LargeCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    x_mw = cloudMW.mw.event.x
    x_cl = cloudMW.cloud.event.x
    v_mw = cloudMW.mw.event.v
    v_cl = cloudMW.cloud.event.v
    u₀ = SVector{12,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3], x_cl[1], x_cl[2], x_cl[3],
                                    v_mw[1], v_mw[2], v_mw[3], v_cl[1], v_cl[2], v_cl[3])
    prob = ODEProblem(ode_perf, u₀, t_span, cloudMW)
    sol  = solve(prob, solver; options...)
    sys_orb = Vector{Orbit}(undef, 2)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[7:9,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[10:12,:])
    return sys_orb
end


# aca estoy...
"""Evolution of a LargeCloudMW (<: GalacticSystem) with dynamical friction"""
function evolve(fric::F, cloudMW::LargeCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F<:AbstractFriction, R<:Real}
    x_mw = cloudMW.mw.event.x
    x_cl = cloudMW.cloud.event.x
    v_mw = cloudMW.mw.event.v
    v_cl = cloudMW.cloud.event.v
    u₀ = SVector{12,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3], x_cl[1], x_cl[2], x_cl[3],
                                    v_mw[1], v_mw[2], v_mw[3], v_cl[1], v_cl[2], v_cl[3])
    p = (fric, cloudMW)
    prob = ODEProblem(ode, u₀, t_span, p)
    sol  = solve(prob, solver; options...)
    sys_orb = Vector{Orbit}(undef, 2)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[7:9,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[10:12,:])
    return sys_orb
end