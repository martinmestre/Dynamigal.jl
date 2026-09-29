"""Potential functions"""


"""
    potential(pot::P, x::Vector{<:Unitful.Length}, t::T) where {P<:AbstractPotential, T<:Unitful.Time}
Unitful Potential for AbstractPotentials con y sin dependencia tempora.
Sirve incluso para AbstractCompositePotential.
"""
function potential(pot::P, x::Vector{<:Unitful.Length}, t::T) where {P<:AbstractPotential, T<:Unitful.Time}
    x, t = adimensional(x, t)
    return potential(pot, x, t)*𝕦.p
end
function potential(pot::P, x::Vector{<:Unitful.Length}) where {P<:AbstractPotential}
    x = adimensional(x)
    return potential(pot, x)*𝕦.p
end

""" Generic method for AbstractCompositePotential
This method is used, for example, by the Exponential3MN potential,
and then in the field <potentials>, contains a CompositePotential.
"""
# potential(pot::C, x::AbstractArray{L}) where {C<:AbstractCompositePotential, L<:Real} =
#     potential(pot.potentials, x)

"""
   potential(pot::C, x::AbstractVector{L}, t::T=0.0) where {C<:AbstractCompositePotential, L<:Real, T<:Real}
Composite Potential
"""
function potential(pot::C, x::AbstractVector{L}, t::T=0.0) where {C<:AbstractCompositePotential, L<:Real, T<:Real}
    sum_pot = zero(L)
    for p ∈ pot
        sum_pot += potential(p, x, t)
    end
    return sum_pot
end


"""
potential(pot::P, x::AbstractVector{L}, t::T) where {P<:AbstractStaticPotential, L<:Real, T<:Real}
Bridge function for static potentials
"""
function potential(pot::P, x::AbstractVector{L}, t::T) where {P<:AbstractStaticPotential, L<:Real, T<:Real}
    return potential(pot, x)
end

"""
potential(pot::P, x::AbstractVector{L}) where {P<:AbstractStaticPotential, L<:Real}
Bridge function for spherical static potentials
"""
function potential(pot::P, x::AbstractVector{L}) where {P<:AbstractSphericalStaticPotential, L<:Real}
    r  = sqrt(  dot(x,x)  )
    return potential(pot, r)
end


""" Potential corresponding to a single spherical and static potential, from ρ(r),
    considering  Φ(∞) = 0.
    Φ(r) = -4πG*[ (1/r) ∫₀ʳ ρ(s) s² ds + ∫ᵣ^∞ ρ(s)s ds ]
    Used only when setting the corresponding specific method:
    potential(pot::PowerLawCutoff, r::L) where {L<:Real} = potential(DensityTrait(), pot, r)
"""
function potential(::DensityTrait, pot::P, r::L) where {P<:AbstractSphericalStaticPotential, L<:Real}
    f(s,p) = s * density(pot, s)
    prob = IntegralProblem(f, (r, Inf))
    return  -G * ( mass(DensityTrait(), pot, r) / r + 4π*solve(prob, QuadGKJL()).u )
end

""" Potential corresponding to a single spherical and static potential, from a(r),
    considering  Φ(∞) = 0.
    Φ(r) = ∫ᵣ^∞ a(s) ds
    Used only when setting the corresponding specific method:
    potential(pot::PowerLawCutoff, r::L) where {L<:Real} = potential(AccelerationTrait(), pot, r)
"""
function potential(::AccelerationTrait, pot::P, r::L) where {P<:AbstractSphericalStaticPotential, L<:Real}
    f(s,p) = acceleration(pot, s)
    prob = IntegralProblem(f, (r, Inf))
    return solve(prob, QuadGKJL()).u
end


"""Concrete potentials"""

"""Spherical"""

"""Allen and Santillan (generalized) halo
Corresponds to Eq. (20) at Irrgang et al. (2013).
The exact A&S corresponds exactly to the case Λ=100kpc and γ=2.02.
"""
function potential(pot::AllenSantillanHalo, r::L) where {L<:Real}
    @unpack_AllenSantillanHalo pot
    f(y) = 1 + (y/a)^(γ-1)
    if r < Λ
        println("r<Λ")
        res = G*(m/a)*( log(f(r)/f(Λ))/(γ-1) - (1-1/f(Λ)) )
    else
        println("else")
        res = -G*(m/r)*(Λ/a)^γ/f(Λ)
    end
    return res
end


"""Hernquist potential"""
function potential(pot::Hernquist, r::L) where {L<:Real}
    @unpack m, a = pot
    return -G*m / (a + r)
end


"""Kepler potential"""
function potential(pot::Kepler, x::AbstractVector{L}) where {L<:Real}
    return -G*pot.m / r
end


"""NFW halo potential"""
function potential(pot::NFW, r::L) where {L<:Real}
    @unpack m, a = pot
    return -G*m*log(1+r/a)/r
end


"""Oscillatory Kepler dependent"""
function potential(pot::OscillatoryKepler, x::AbstractVector{L}, t::T) where {L<:Real, T<:Real}
    @unpack m, τ = pot
    return -G*m*sin((2π/τ)*t) / sqrt(t^2 +  dot(x,x) )
end


"""Plummer potential"""
function potential(pot::Plummer, r::L) where {L<:Real}
     @unpack m, a = pot
    return -G*m / sqrt(a^2 +  r^2)
end

"""PowerLawCutoff potential -- to do: compute analytically"""
potential(pot::PowerLawCutoff, r::L) where {L<:Real} = potential(AccelerationTrait(), pot, r)


"""Axisymmetric"""

"""Miyamoto-Nagai disk potential"""
function potential(pot::MiyamotoNagai, x::AbstractArray{L}) where {L<:Real}
    @unpack m, a, b = pot
    y = @view x[1:2]
    return -G*m/sqrt( dot(y,y) + (a + sqrt(b^2+x[3]^2))^2 )
end

