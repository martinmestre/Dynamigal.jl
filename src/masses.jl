
""" Enclosed mass for P <:AbstractSphericalStaticPotential at position x"""
function mass(pot::P, x::AbstractVector{L}) where {P<:AbstractSphericalStaticPotential, L<:Real}
    r = sqrt( dot(x,x) )
    return mass(pot, r)
end


"""NFW mass"""
mass(pot::NFW, r::L) where {L<:Real} = pot.m*f_nfw(r/pot.a)

"""
    PowerLawCutoff mass(r)

    density(pot::PowerLawCutoff, r::AbstractVector{L}) where {L<:Real}
    Expression from Gala.
"""
function mass(pot::PowerLawCutoff, r::L) where {L<:Real}
    @unpack_PowerLawCutoff pot
    return 2π*𝔸*c^(3-α)*gamma(β)*gamma_inc(β, r*r/(c*c),0)[1]
end