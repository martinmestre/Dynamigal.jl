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
        return 0.0
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
        return 0.0
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
        return 0.0
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
        return 0.0
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
        return 0.0
    else
        Λ = r / γₕ / max(rₚ, G*mₚ/ν²)
        lnΛ = 0.5*log(1+Λ^2)
        χ = ν /(σₕ(r)*√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end
