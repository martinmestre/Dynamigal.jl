"""Potential functions"""

"""General trait dispatch"""
function potential(pot::P, x::AbstractArray{T}, t::D) where {P<:AbstractPotential, T, D}
    return potential(time_dependence(P), pot, x, t)
end

function potential(::TimeIndependent, pot::P, x::AbstractArray{T}, t::D) where {P<:AbstractPotential, T, D}
    return potential(pot, x)
end
function potential(::TimeDependent, pot::P, x::AbstractArray{T}, t::D) where {P<:AbstractPotential, T, D}
    return potential(pot, x, t)
end

"""Unitful Potential for AbstractPotentials"""
function potential(pot::P, x::Vector{<:Unitful.Length}, t::T) where {P<:AbstractPotential, T<:Unitful.Time}
    x, t = adimensional(x, t)
    return potential(pot, x, t)*𝕦.p
end
function potential(pot::P, x::Vector{<:Unitful.Length}) where {P<:AbstractPotential}
    x = adimensional(x)
    return potential(pot, x)*𝕦.p
end

"""Composite Potential"""
function potential(pot::CompositePotential, x::AbstractArray{L}, t::T) where {L<:Real, T<:Real}
    sum_pot = zero(L)
    for p ∈ pot
        sum_pot += potential(p, x, t)
    end
    return sum_pot
end


"""List of specific Potentials..."""


"""Kepler potential"""
function potential(pot::Kepler, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / sqrt(x'x)
end

"""Kepler time dependent"""
function potential(pot::OscillatoryKepler, x::AbstractArray{L}, t::T) where {L<:Real, T<:Real}
    @unpack m, τ = pot
    return -G*m*sin((2π/τ)*t) / sqrt(t^2 + x'x)
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
