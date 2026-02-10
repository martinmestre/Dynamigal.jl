"""Dissipative Forces"""


"""Chandrasekhar dynamical friction configuration"""
@with_kw struct ChandrasekharFriction{T<:Real, R<:Real, S<:Real} <:AbstractFriction
    mₚ::R # perturber mass
    lnΛ::T  # Coulomb logarithm
    σₕ::S # host's mean velocity dispersion
end
ChandrasekharFriction(lnΛ::T, mₚ::R, σₕ::S) where {T<:Real, R<:Unitful.Mass, S<:Unitful.Velocity} =
    ChandrasekharFriction(lnΛ, ustrip(uconvert(𝕦.m, mₚ)),  ustrip(uconvert(𝕦.v, σₕ)) )

"""Galpy's Chandrasekhar dynamical friction configuration"""
@with_kw struct GalpyFriction{M<:Real, L<:Real, V<:Real, R<:Real} <:AbstractFriction
    mₚ::M # perturber mass
    rₚ::L # perturber's half-mass radius.
    σₕ::V # host's mean velocity dispersion
    γₕ::R = 1.0  # according to DOI 10.1093/mnras/stw2011, this quantity should be
                # γ = max(|r/ρ * dρ/dr|, r)
                # I demonstrated that the equation in Galpy manual (https://docs.galpy.org/en/v1.11.0/ reference/potentialchandrasekhardynfric.html#dynamfric-potential)
                # is equivalent to Eq. (6) in reference DOI 10.1093/mnras/stw2011.
end

"""Agama's Chandrasekhar dynamical friction configuration"""
@with_kw struct AgamaFriction{M<:Real, L<:Real, V<:Real} <:AbstractFriction
    mₚ::M # perturber mass
    rₚ::L = 0.85*(mₚ/1.0e11)^0.6 # perturber's scale radius
    bₘ::L = 2rₚ # host's mean velocity dispersion
    σₕ::V # host's mean velocity dispersion
end