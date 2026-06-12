"""Build friction arrays"""
function build_friction(mps::MacroParticleSystem, algorithm::F) where {F}
    @warn "The value of α is not exact for general potentials"
    n = length(mps)
    α = 1+sqrt(2)
    γ = 0.1
    η = 100.0
    fric = Matrix{algorithm}(undef, n, n)
    for j in 1:n
        if mps[j].pot isa CompositePotential
            r_min = γ*minimum([mps[j].pot[k].a for k in 1:length(mps[j].pot)])
            r_max = η*maximum([mps[j].pot[k].a for k in 1:length(mps[j].pot)])
        else
            r_min = γ*mps[j].pot.a
            r_max = η*mps[j].pot.a
        end
        σ=velocity_dispersion(mps[j].pot; r_min=r_min, r_max=r_max)
        for i in 1:n
            if mps[i].pot isa CompositePotential
                mass = sum(mps[i].pot[k].m for k in 1:length(mps[i].pot))
                scale = maximum([mps[i].pot[k].a for k in 1:length(mps[i].pot)])
            else
                mass = mps[i].pot.m
                scale = mps[i].pot.a
            end
            fric[i,j] = algorithm(mass, α*scale, σ)
        end
    end
    return fric
end

function build_friction!(fric::Matrix{GalpyFriction}, mps::MacroParticleSystem)
    @warn "The value of α is not exact for general potentials"
    n = length(mps)
    α = 1+sqrt(2)
    γ = 0.1
    η = 100.0
    for j in 1:n
        if mps[j].pot isa CompositePotential
            r_min = γ*minimum([mps[j].pot[k].a for k in 1:length(mps[j].pot)])
            r_max = η*maximum([mps[j].pot[k].a for k in 1:length(mps[j].pot)])
        else
            r_min = γ*mps[j].pot.a
            r_max = η*mps[j].pot.a
        end
        σ=velocity_dispersion(mps[j].pot; r_min=r_min, r_max=r_max)
        for i in 1:n
            if mps[i].pot isa CompositePotential
                mass = sum(mps[i].pot[k].m for k in 1:length(mps[i].pot))
                scale = maximum([mps[i].pot[k].a for k in 1:length(mps[i].pot)])
            else
                mass = mps[i].pot.m
                scale = mps[i].pot.a
            end
            fric[i,j] = GalpyFriction(mass, α*scale, σ)
        end
    end
end

"""Build friction arrays"""
function build_friction_pyramid!(fric::Matrix{GalpyFriction}, mps::MacroParticleSystem)
    @warn "The value of α is not exact for general potentials. Partial friction.."
    n = length(mps)
    α = 1+sqrt(2)
    γ = 0.1
    η = 100.0
    for j in 1:n
        if mps[j].pot isa CompositePotential
            r_min = γ*minimum([mps[j].pot[k].a for k in 1:length(mps[j].pot)])
            r_max = η*maximum([mps[j].pot[k].a for k in 1:length(mps[j].pot)])
            mass_host = sum(mps[j].pot[k].m for k in 1:length(mps[j].pot))
        else
            r_min = γ*mps[j].pot.a
            r_max = η*mps[j].pot.a
            mass_host = mps[j].pot.m
        end
        σ=velocity_dispersion(mps[j].pot; r_min=r_min, r_max=r_max)
        for i in 1:n
            if mps[i].pot isa CompositePotential
                mass = sum(mps[i].pot[k].m for k in 1:length(mps[i].pot))
                scale = maximum([mps[i].pot[k].a for k in 1:length(mps[i].pot)])
            else
                mass = mps[i].pot.m
                scale = mps[i].pot.a
            end
            if mass_host >= mass
                fric[i,j] = GalpyFriction(mass, α*scale, σ)
            else
                fric[i,j] = GalpyFriction(0.0, α*scale, σ)
            end
        end
    end
end

function build_friction(host::T, satellite::S, algorithm::F) where {T<:AbstractMacroParticle,S<:AbstractMacroParticle, F}
    @warn "The value of α is not exact for general potentials"
    α = 1+sqrt(2)
    γ = 0.1
    η = 50.0
    if host.pot isa CompositePotential
        r_min = γ*minimum([host.pot[k].a for k in 1:length(host.pot)])
        r_max = η*maximum([host.pot[k].a for k in 1:length(host.pot)])
    else
        r_min = γ*host.pot.a
        r_max = η*host.pot.a
    end

    if satellite.pot isa CompositePotential
        mass = sum(satellite.pot[k].m for k in 1:length(satellite.pot))
        scale = maximum([satellite.pot[k].a for k in 1:length(satellite.pot)])
    else
        mass = satellite.pot.m
        scale = satellite.pot.a
    end
    # @show r_min r_max scale mass algorithm
    σ=velocity_dispersion(host.pot; r_min=r_min, r_max=r_max)

    if algorithm == GalpyFriction
        return algorithm(mass, α*scale, σ)
    elseif algorithm == TangoFriction
        return algorithm(mass, σ)
    elseif algorithm == AgamaFriction
        return algorithm(mass, σ)
    else
        @error "Friction algorithm incorrect"
    end
end


