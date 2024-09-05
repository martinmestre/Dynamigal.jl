"""Potential functions"""


"""Unitful Potential for UnionAbstractPotentials"""
function potential(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, t::T) where {T<:Unitful.Time}
    x, t = adimensional(x, t)
    return potential(pot, x, t)*𝕦.p
end


"""Potential of a sum of AbstractPotentials"""
function potential(pot::Vector{<:AbstractPotential}, x::AbstractArray{L}, t::T) where {L<:Real, T<:Real}
    sum_pot = zeros(3)
    for i ∈ eachindex(pot)
        sum_pot .+= potential(pot[i], x, t)
    end
    return sum_pot
end

"""Multiple dispatch when t variable is/is-not used"""
potential(pot::UnionAbstractPotentials, x::AbstractArray{T}, t::D) where {T<:Real, D<:Real} =
    potential(pot, x)
potential(pot::UnionAbstractPotentials, x::AbstractArray{<:Unitful.Length}) =
    potential(pot, x, 0𝕦.t)


    """List of specific Potentials..."""

"""TimeDependent potential"""
function potential(pot::TimeDependent, x::AbstractArray{L}, t::T) where {L<:Real, T<:Real}
    return -G*pot.m / sqrt(t^2+x'x)
end

"""Kepler potential"""
function potential(pot::Kepler, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / sqrt(x'x)
end

"""Plummer potential"""
function potential(pot::Plummer, x::AbstractArray{T}) where {T<:Real}
    @unpack m, a = pot
    return -G*m / sqrt(a^2 + x'x)
end

"""Hernquist potential"""
function potential(pot::Hernquist, x::AbstractArray{T}) where {T<:Real}
    @unpack m, a = pot
    return -G*m / (a + sqrt(x'x))
end

"""Miyamoto-Nagai disk potential"""
function potential(pot::MiyamotoNagaiDisk, x::AbstractArray{T}) where {T<:Real}
    @unpack m, a, b = pot
    return -G*m/sqrt( x[1:2]'x[1:2] + (a + sqrt(b^2+x[3]^2))^2 )
end

"""Allen and Santillan (generalized) halo"""
function potential(pot::AllenSantillanHalo, x::AbstractArray{T}) where {T<:Real}
    @unpack_AllenSantillanHalo pot
    f(y) = 1 + (y/a)^(γ-1)
    r  = sqrt( x'x )
    if r < Λ
        res = -G*(m/a)*( log(f(r)/f(Λ))/(γ-1) - (1-1/f(Λ)) )
    else
        res = -G*(m/r)*(Λ/a)^γ/f(Λ)
    end
    return res
end

"""NFW halo potential"""
function potential(pot::NFW, x::AbstractArray{T}) where {T<:Real}
    @unpack m, a, 𝔸 = pot
    r = sqrt(x'x)
    return -G*m/𝔸*log(1+r/a)/r
end
