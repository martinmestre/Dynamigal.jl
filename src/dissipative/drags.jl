"""
Dynamical friction acceleration for standar Chandrasekhar formula.
x = position relative to potencial source centre.
v = velocity relative to potential source velocity.
"""
function drag(fric::ChandrasekharFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack lnΛ, mₚ, σₕ = fric
    ν = sqrt(dot(v,v))
    if ν < 𝕗.ϵ_ν
        return 0.0
    else
        χ = ν /(σₕ√2)
        return -4π*G*lnΛ*density(p, x, t) * mₚ * ( erf(χ) - (2/√π)*χ*exp(-χ^2) ) * v/ν^3
    end
end