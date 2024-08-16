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
    return -G*pot.m / sqrt(pot.a^2 + x'x)
end

"""Hernquist potential"""
function potential(pot::Hernquist, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / (pot.a + sqrt(x'x))
end

"""Miyamoto-Nagai disk potential"""
function potential(pot::MiyamotoNagaiDisk, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m/sqrt( x[1:2]'x[1:2] + (pot.a + sqrt(pot.b^2+x[3]^2))^2 )
end

"""Allen and Santillan (generalized) halo"""
function potential(pot::AllenSantillanHalo, x::AbstractArray{T}) where {T<:Real}
    f(y) = 1 + (y/pot.a)^(pot.γ-1)
    r  = sqrt( x'x )
    if r < pot.Λ
        res = -G*(pot.m/pot.a)*( log(f(r)/f(pot.Λ))/(pot.γ-1) - (1-1/f(pot.Λ)) )
    else
        res = -G*(pot.m/r)*(pot.Λ/pot.a)^pot.γ/f(pot.Λ)
    end
    return res
end

"""NFW halo potential"""
@inline f_nfw(x) = log(1+x)-x/(1+x)

function potential(pot::NFW, x::AbstractArray{T}) where {T<:Real}
    𝔸 = f_nfw(concentration(pot))
    r = sqrt(x'x)
    return -G*pot.m/𝔸*log(1+r/pot.a)/r
end
