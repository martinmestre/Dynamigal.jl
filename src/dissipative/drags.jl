"""
Dynamical friction acceleration for standar Chandrasekhar formula.
x = position relative to potencial source centre.
v = velocity relative to potential source velocity.
"""
function drag(fric::ChandrasekharFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, lnΛ, σₕ = fric
    ν = sqrt(dot(v,v))
    if ν < 𝕗.ϵ_ν
        return 0.0
    else
        χ = ν /(σₕ√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end

function drag(fric::GalpyFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, rₚ, σₕ, γₕ = fric
    r = sqrt(dot(x,x))

    ν² = dot(v,v)
    ν = sqrt(ν²)

    if ν < 𝕗.ϵ_ν
        return 0.0
    else
        Λ = r / γₕ / max(rₚ, G*mₚ/ν²)
        lnΛ = 0.5*log(1+Λ^2)
        χ = ν /(σₕ√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end

function drag(fric::AgamaFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack mₚ, rₚ, bₘ, σₕ = fric
    r = sqrt(dot(x,x))
    ν = sqrt(dot(v,v))
    lnΛ = max(0, sqrt( log(r/bₘ) ) )
    if ν < 𝕗.ϵ_ν
        return 0.0
    else
        χ = ν /(σₕ√2)
        return -4π*G^2*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end