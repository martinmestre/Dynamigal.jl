"""
Dynamical friction acceleration for standard Chandrasekhar formula.
x = position relative to potencial source centre.
v = velocity relative to potential source velocity.
"""




"""Chandrasekhar dynamical friction drag
Exact formula from Binney & Tremaine (2008). No assumptions on b_max or v_typ.
"""
function drag(fric::ChandrasekharFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, rₚ, b_max, v_typ, σₕ = fric
    ν = sqrt(dot(v,v))
    if ν < 𝕗.ϵ_ν
        return zeros(3)
    else
        lnΛ = log( b_max / max(rₚ, G*mₚ/v_typ^2) )
        χ = ν /(σₕ√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end

"""Gala's Chandrasekhar dynamical friction drag
Identical to Adrian's hack for Gala:
https://gist.github.com/adrn/2ca7ed34b38a47252dfbecdb0c70bd97
∃ assumptions on b_max=|x| and v_typ=σₕ.
"""
function drag(fric::GalaFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, rₚ, σₕ = fric
    r = sqrt(dot(x,x))
    ν = sqrt(dot(v,v))
    if ν < 𝕗.ϵ_ν
        return zeros(3)
    else
        lnΛ = log( r / max(rₚ, G*mₚ/σₕ(r)^2) )
        χ = ν /(σₕ(r)*√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end


"""Tango for three's Chandrasekhar dynamical friction drag
Taken from MNRAS 501, 2279–2304 (2021), Vasiliev et al., page 2285.
The lnΛ=constant recipe is only used for Sagittarius dwarf, not for the clouds, so
this is not similar to Agama's formula below.
Besides, σ=constant.
"""
function drag(fric::TangoFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, lnΛ, σₕ = fric
    ν = sqrt(dot(v,v))
    if ν < 𝕗.ϵ_ν
        return zeros(3)
    else
        χ = ν /(σₕ√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end


"""Agama's Chandrasekhar dynamical friction configuration
See formula here:
https://github.com/GalacticDynamics-Oxford/Agama/blob/1a0c519f3d89c621f04c2e4502183e22dc7e441a/py/example_lmc_mw_interaction.py#L96
👀 Note the √ in the lnΛ formula!
"""
function drag(fric::AgamaFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, b_min, σₕ = fric
    r = sqrt(dot(x,x))
    ν = sqrt(dot(v,v))
    lnΛ = max(0, sqrt(log(r / b_min)) )
    if ν < 𝕗.ϵ_ν
        return zeros(3)
    else
        χ = ν /(σₕ(r)*√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end


"""Galpy's Chandrasekhar dynamical friction configuration
See formula here:
https://docs.galpy.org/en/latest/reference/potentialchandrasekhardynfric.html#galpy.potential.ChandrasekharDynamicalFrictionForce
Src code here:
https://github.com/jobovy/galpy/blob/main/galpy/potential/ChandrasekharDynamicalFrictionForce.py#L37-L159
Here we do not have cut the drag when r<r_min or r>r_max, like in Galpy.
"""
function drag(fric::GalpyFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, rₚ, σₕ, γₕ = fric
    r = sqrt(dot(x,x))
    ν² = dot(v,v)
    ν = sqrt(ν²)
    if ν < 𝕗.ϵ_ν
        return zeros(3)
    else
        Λ = r / γₕ / max(rₚ, G*mₚ/ν²)
        lnΛ = 0.5*log(1+Λ^2)
        χ = ν /(σₕ(r)*√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end

"""Build friction arrays"""
function build_friction!(fric::Matrix{F}, mps::MacroParticleSystem) where {F<:AbstractFriction}
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
            fric[i,j] = F(mass, α*scale, σ)
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
    @show r_min r_max scale mass algorithm
    σ=velocity_dispersion(host.pot; r_min=r_min, r_max=r_max)

    return algorithm(mass, α*scale, σ)
end


