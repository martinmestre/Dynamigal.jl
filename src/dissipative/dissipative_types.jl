"""Dissipative Forces"""


"""Chandrasekhar dynamical friction configuration"""
@with_kw struct ChandrasekharFriction{T<:Real, R<:Real, S<:Real} <:AbstractFriction
    lnΛ::T  # Coulomb logarithm
    mₚ::R # perturber mass
    σₕ::S # host's mean velocity dispersion
end
ChandrasekharFriction(lnΛ::T, mₚ::R, σₕ) where {T<:Real, R<:Unitful.Mass, S<:Unitful.Velocity} =
    ChandrasekharFriction(lnΛ, ustrip(uconvert(𝕦.m, mₚ)),  ustrip(uconvert(𝕦.v, σₕ)) )