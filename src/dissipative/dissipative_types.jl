"""Dissipative Forces"""



"""Chandrasekhar dynamical friction configuration
Exact formula from Binney & Tremaine (2008). No assumptions on b_max or v_typ.
"""
@with_kw struct ChandrasekharFriction{M<:Real, L<:Real, D<:Real, V<:Real, R<:Real} <:AbstractFriction
    mₚ::M # perturber mass
    rₚ::L # perturbers half-mass radius (=r_hm in B&T 2008)
    b_max::D  # maximum orthogonal impact distance
    v_typ::V # typical speed in the host
    σₕ::R # host's mean velocity dispersion
end

ChandrasekharFriction(mₚ::M, rₚ::L, b_max::D, v_typ::V, σₕ::R) where {M<:Unitful.Mass, L<:Unitful.Length,  D<:Unitful.Length, V<:Unitful.Velocity, R<:Unitful.Velocity} =
    ChandrasekharFriction(ustrip(uconvert(𝕦.m, mₚ)), ustrip(uconvert(𝕦.l, rₚ)),
    ustrip(uconvert(𝕦.l, b_max)), ustrip(uconvert(𝕦.v, v_typ)), ustrip(uconvert(𝕦.v, σₕ)) )


"""Gala's Chandrasekhar dynamical friction configuration
Simular to Adrian's hack for Gala:
https://gist.github.com/adrn/2ca7ed34b38a47252dfbecdb0c70bd97
∃ assumptions on b_max=|x| and v_typ=σₕ,
but using a callable σ parameter: σ(r).
"""
@with_kw struct GalaFriction{M<:Real, L<:Real, V} <:AbstractFriction
    mₚ::M # perturber mass
    rₚ::L  # perturbers half-mass radius (=r_hm in B&T 2008)
    σₕ::V # host's mean velocity dispersion
end

GalaFriction(mₚ::M, rₚ::L, σₕ::V) where {M<:Unitful.Mass, L<:Unitful.Length, V} =
    GalaFriction(ustrip(uconvert(𝕦.m, mₚ)), ustrip(uconvert(𝕦.l, rₚ)), σₕ )


"""Tango for three's Chandrasekhar dynamical friction configuration
Taken from MNRAS 501, 2279–2304 (2021), Vasiliev et al., page 2285.
The lnΛ=constant recipe is only used for Sagittarius dwarf, not for the clouds, so
this is not similar to Agama's formula below.
Besides, σ=constant.
"""
@with_kw struct TangoFriction{T<:Real, R<:Real, S<:Real} <:AbstractFriction
    mₚ::R # perturber mass
    lnΛ::T  # Coulomb logarithm
    σₕ::S # host's mean velocity dispersion
end
TangoFriction(mₚ::R, lnΛ::T, σₕ::S) where {R<:Unitful.Mass, T<:Real, S<:Unitful.Velocity} =
    TangoFriction(ustrip(uconvert(𝕦.m, mₚ)), lnΛ, ustrip(uconvert(𝕦.v, σₕ)) )


"""Agama's Chandrasekhar dynamical friction configuration
See formula here:
https://github.com/GalacticDynamics-Oxford/Agama/blob/1a0c519f3d89c621f04c2e4502183e22dc7e441a/py/example_lmc_mw_interaction.py#L96
"""
struct AgamaFriction{M<:Real, L<:Real, T<:Real, F} <:AbstractFriction
    mₚ::M # perturber mass
    rₚ::L  # perturber's radius scale
    b_min::T # minimum impact parameter in the Coulomb logarithm
    σₕ::F # host's mean velocity dispersion
end
function AgamaFriction(mₚ::M, σₕ::F) where {M<:Real, F}  # This is the prescription in Agama script.
    rₚ = 0.85*(mₚ/1.0e11)^0.6
    b_min = 2rₚ
    return AgamaFriction(mₚ, rₚ, b_min, σₕ)
end
AgamaFriction(mₚ::M, σₕ::F) where {M<:Unitful.Mass, F} =
    AgamaFriction(ustrip(uconvert(𝕦.m, mₚ)), σₕ)


"""Galpy's Chandrasekhar dynamical friction configuration
See formula here:
https://docs.galpy.org/en/latest/reference/potentialchandrasekhardynfric.html#galpy.potential.ChandrasekharDynamicalFrictionForce
Src code here:
https://github.com/jobovy/galpy/blob/main/galpy/potential/ChandrasekharDynamicalFrictionForce.py#L37-L159
"""
struct GalpyFriction{M<:Real, L<:Real, F, R<:Real} <:AbstractFriction
    mₚ::M # perturber mass
    rₚ::L # perturber's half-mass radius.
    σₕ::F # host's mean velocity dispersion
    γₕ::R # according to DOI 10.1093/mnras/stw2011, this quantity should be
                # γ = max(|r/ρ * dρ/dr|, r)
                # I demonstrated that the equation in Galpy manual (https://docs.galpy.org/en/v1.11.0/ reference/potentialchandrasekhardynfric.html#dynamfric-potential)
                # is equivalent to Eq. (6) in reference DOI 10.1093/mnras/stw2011.
end
function GalpyFriction(mₚ::M, rₚ::L, σₕ::F) where {M<:Real, L<:Real, F}
    return GalpyFriction(mₚ, rₚ, σₕ, 1.0)
end
GalpyFriction(mₚ::M, rₚ::L, σₕ::F) where {M<:Unitful.Mass, L<:Unitful.Length, F} =
    GalpyFriction(ustrip(uconvert(𝕦.m, mₚ)), ustrip(uconvert(𝕦.l, rₚ)), σₕ)


