"""
Dynamical friction acceleration for standar Chandrasekhar formula.
x = position relative to potencial source centre.
v = velocity relative to potential source velocity.
"""
function drag(𝕗::ChandrasekharFriction, p::P, x::AbstractArray{L}, v::AbstractArray{L}, t::T) where {P<:AbstractPotential, L<:Real, T<:Real}
    @unpack lnΛ, mₚ, σₕ = 𝕗
    ν = sqrt(dot(v,v))
    χ = ν /(σₕ√2)
    return -4π*G*lnΛ*density(p, x, t) * mₚ * v/ν^3 * ( erf(χ) - (2/√π)*χ*exp(-χ^2) )
end