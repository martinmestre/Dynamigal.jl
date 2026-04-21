"""Evolution functions"""


"""Evolution of a an initial condition in an AbstractPotential"""
function _evolve(pot::P, x::AbstractVector{D}, v::AbstractVector{F},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential, D<:Real, F<:Real, T}
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(_ode, u₀, t_span, p)
    sol = solve(prob, solver; options...)
    orb = Orbit(sol.t, sol[sis,:], sol[siss,:])
    return orb
end
# evolve(below) is a little faster than _evolve(above)
# ojo! en test de Orbits.ipynb da mejor la de arriba
function evolve(pot::P, x::AbstractVector{D}, v::AbstractVector{F},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {P<:AbstractPotential, D<:Real, F<:Real, T}
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
    sol = solve(prob, solver; options...)
    orb = Orbit(sol.t, sol[sis,:], sol[siss,:])
    return orb
end

"""Evolution of a an initial condition in an AbstractPotential including dynamical friction"""
function evolve(fric::F, pot::P, x::AbstractVector{L}, v::AbstractVector{V},
   t_span::Tuple{T,T}, solver=𝕤.ode; options=ntSolverOptions()) where {F<:AbstractFriction, P<:AbstractPotential, L, V, T}
    p = (fric, pot)
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
function evolve(::RawSolutionTrait, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, mps)
    return solve(prob, solver; options...)
end

"""Evolution of a system of MacroParticle, general type"""
function evolve(::GenSysTrait, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    println("🎲 Evolving a MPS system without friction, GenSystTrait")
    sol  = evolve(RawSolutionTrait(), mps, t_span, solver; options=options)
    sys_orb = Vector{Orbit}(undef, length(mps))
    n = length(sol.u[begin])÷2
    for i ∈ eachindex(mps)
        j_x = selec(i)
        j_v = n+j_x
        sys_orb[i] = Orbit(sol.t, sol[j_x:j_x+2,:], sol[j_v:j_v+2,:])
    end
    return sys_orb
end


"""Evolution of a system of MacroParticles using Mutating ODE trait"""
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

"""Evolution of a system of MacroParticle with mutual dynamical friction considered.
The integration scheme is similar to ::GenSysTrait above but with the friction force."""
function evolve(::MutualFrictionTrait, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
    u₀ = SA[x...,v...]
    n = length(mps)
    fric = Matrix{GalpyFriction}(undef, n, n)
    build_friction!(fric, mps)
    p = (fric, mps)
    prob = ODEProblem(ode, u₀, t_span, p)
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

"""Evolution of a system of MacroParticle with partial (or full depending on fric matrix)
dynamical friction considered. The integration scheme is similar to ::GenSysTrait above but with the friction force."""
function evolve(::RawSolutionTrait, fric::Matrix{F}, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F<:AbstractFriction, R<:Real}
    x = reduce(vcat, [mps[i].event.x for i in eachindex(mps)])
    v = reduce(vcat, [mps[i].event.v for i in eachindex(mps)])
    u₀ = SA[x...,v...]
    p = (fric, mps)
    prob = ODEProblem(ode, u₀, t_span, p)
    return solve(prob, solver; options...)
end
function evolve(fric::Matrix{F}, mps::MacroParticleSystem, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F<:AbstractFriction, R<:Real}
    sol  = evolve(RawSolutionTrait(), fric, mps, t_span, solver; options=options)
    sys_orb = Vector{Orbit}(undef, length(mps))
    n = length(sol.u[begin])÷2
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


"""Evolution of a LargeCloudMW (<: GalacticSystem) with dynamical friction for the cloud"""
function evolve(::RawSolutionTrait, fric::F, cloudMW::LargeCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F<:AbstractFriction, R<:Real}
    x_mw = cloudMW.mw.event.x
    x_cl = cloudMW.cloud.event.x
    v_mw = cloudMW.mw.event.v
    v_cl = cloudMW.cloud.event.v
    u₀ = SVector{12,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3], x_cl[1], x_cl[2], x_cl[3],
                                    v_mw[1], v_mw[2], v_mw[3], v_cl[1], v_cl[2], v_cl[3])
    p = (fric, cloudMW)
    prob = ODEProblem(ode, u₀, t_span, p)
    return solve(prob, solver; options...)
end
function evolve(fric::F, cloudMW::LargeCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F<:AbstractFriction, R<:Real}
    sol  = evolve(RawSolutionTrait(), fric, cloudMW, t_span, solver; options=options)
    sys_orb = Vector{Orbit}(undef, 2)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[7:9,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[10:12,:])
    return sys_orb
end

"""Evolution of a CloudsMW (<: GalacticSystem) without friction."""
function evolve(::RawSolutionTrait, system::CloudsMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    # println("🎲 Evolving a CloudsMW system without friction, RawSolutionTrait")
    x_mw = system.mw.event.x
    x_lc = system.large.event.x
    x_sc = system.small.event.x
    v_mw = system.mw.event.v
    v_lc = system.large.event.v
    v_sc = system.small.event.v
    u₀ = SVector{18,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3],
                                    x_lc[1], x_lc[2], x_lc[3],
                                    x_sc[1], x_sc[2], x_sc[3],
                                    v_mw[1], v_mw[2], v_mw[3],
                                    v_lc[1], v_lc[2], v_lc[3],
                                    v_sc[1], v_sc[2], v_sc[3])
    p = system
    prob = ODEProblem(ode, u₀, t_span, p)
    return solve(prob, solver; options...)
end
function evolve(system::CloudsMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {R<:Real}
    sol = evolve(RawSolutionTrait(), system, t_span, solver; options=options)
    sys_orb = Vector{Orbit}(undef, 3)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[10:12,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[13:15,:])
    sys_orb[3] = Orbit(sol.t, sol[7:9,:], sol[16:18,:])
    return sys_orb
end

"""Evolution of a CloudsMW (<: GalacticSystem) system including MW's dynamical friction on both clouds and LC's dynamical friction on the SC. And the reflex acceleration of both clouds on the MW."""
function evolve(::RawSolutionTrait, fric::F, system::CloudsMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F, R<:Real}
    # println("🎲 Evolving a CloudsMW system with pyramidal frictions, RawSolutionTrait")
    x_mw = system.mw.event.x
    x_cl = system.large.event.x
    x_sat = system.small.event.x
    v_mw = system.mw.event.v
    v_cl = system.large.event.v
    v_sat = system.small.event.v
    u₀ = SVector{18,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3],
                                    x_cl[1], x_cl[2], x_cl[3],
                                    x_sat[1], x_sat[2], x_sat[3],
                                    v_mw[1], v_mw[2], v_mw[3],
                                    v_cl[1], v_cl[2], v_cl[3],
                                    v_sat[1], v_sat[2], v_sat[3])
    p = (fric, system)
    prob = ODEProblem(ode, u₀, t_span, p)
    sol  = solve(prob, solver; options...)
    return sol
end
function evolve(fric::F, system::CloudsMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F, R<:Real}
    # println("🎲 Evolving a CloudsMW system with pyramidal frictions")
    sol  = evolve(RawSolutionTrait(), fric, system, t_span, solver; options=options)
    sys_orb = Vector{Orbit}(undef, 3)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[10:12,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[13:15,:])
    sys_orb[3] = Orbit(sol.t, sol[7:9,:], sol[16:18,:])
    return sys_orb
end

"""Evolution of a SatelliteCloudMW (<: GalacticSystem) with dynamical friction for the cloud. Two alterantives of friction treatment in acceleration method, depending on the parameter type."""
function evolve(::RawSolutionTrait, fric::F, system::SatelliteCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F, R<:Real}
    # println("🎲 Evolving a SatelliteCloudMW system, RawSolutionTrait")
    x_mw = system.mw.event.x
    x_cl = system.cloud.event.x
    x_sat = system.satellite.event.x
    v_mw = system.mw.event.v
    v_cl = system.cloud.event.v
    v_sat = system.satellite.event.v
    u₀ = SVector{18,typeof(x_mw[1])}(x_mw[1], x_mw[2], x_mw[3],
                                    x_cl[1], x_cl[2], x_cl[3],
                                    x_sat[1], x_sat[2], x_sat[3],
                                    v_mw[1], v_mw[2], v_mw[3],
                                    v_cl[1], v_cl[2], v_cl[3],
                                    v_sat[1], v_sat[2], v_sat[3])
    p = (fric, system)
    prob = ODEProblem(ode, u₀, t_span, p)
    return solve(prob, solver; options...)
end
function evolve(fric::F, system::SatelliteCloudMW, t_span::Tuple{R,R}, solver=𝕤.ode; options=ntSolverOptions()) where {F, R<:Real}
    # println("🎲 Evolving a SatelliteCloudMW system")
    sol  = evolve(RawSolutionTrait(), fric, system, t_span, solver; options=options)
    sys_orb = Vector{Orbit}(undef, 3)
    sys_orb[1] = Orbit(sol.t, sol[1:3,:], sol[10:12,:])
    sys_orb[2] = Orbit(sol.t, sol[4:6,:], sol[13:15,:])
    sys_orb[3] = Orbit(sol.t, sol[7:9,:], sol[16:18,:])
    return sys_orb
end