"""Potential functions"""


"""Unitful Potential for UnionAbstractPotentials"""
function potential(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, t::T=nothing) where {T<:Union{Unitful.Time,Nothing}}
    x, t = code_units(x, t)
    return potential(pot, x, t)*𝕦.p
end


"""Potential of a sum of AbstractPotentials"""
function potential(pot::Vector{<:AbstractPotential}, x::AbstractArray{L}, t::T=nothing) where {L<:Real, T<:Union{Real,Nothing}}
    sum_pot = zeros(3)
    for i ∈ eachindex(pot)
        sum_pot .+= potential(pot[i], x, t)
    end
    return sum_pot
end


"""List of specific Potentials..."""

"""TimeDependent potential"""
function potential(pot::TimeDependent, x::AbstractArray{L}, t::T) where {L<:Real, T<:Real}
    return -G*pot.m / sqrt(t^2+x'x)
end


"""PointMass potential"""
function potential(pot::PointMass, x::AbstractArray{T}, t=nothing) where {T<:Real}
    return -G*pot.m / sqrt(x'x)
end


"""Plummer potential"""
function potential(pot::Plummer, x::AbstractArray{T}, t=nothing) where {T<:Real}
    return -G*pot.m / sqrt(pot.b^2 + x'x)
end


"""Miyamoto-Nagai disk potential"""
function potential(pot::MiyamotoNagaiDisk, x::AbstractArray{T}, t=nothing) where {T<:Real}
    return -G*pot.m/sqrt( x[1:2]'x[1:2] + (pot.a + sqrt(pot.b^2+x[3]^2))^2 )
end


"""Allen and Santillan (generalized) halo"""
function potential(pot::AllenSantillanHalo, x::AbstractArray{T}, t=nothing) where {T<:Real}
    f(y) = 1.0 + (y/pot.a)^(pot.γ-1.0)
    r  = sqrt( x'x )
    if r < pot.Λ
        res = -G*(pot.m/pot.a)*( log(f(r)/f(pot.Λ))/(pot.γ-1.0) - (1.0-1.0/f(pot.Λ)) )
    else
        res = -G*(pot.m/r)*(pot.Λ/pot.a)^pot.γ/f(pot.Λ)
    end
    return res
end

