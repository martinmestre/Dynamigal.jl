"""Potential functions"""


"""Unitful UnionAbstractPotentials"""
function potential(pot::UnionAbstractPotentials, x::Vector{<:U.Length})
    x = ustrip(uconvert.(u_L, x))
    return u_V^2*potential(pot,x)
end


"""Potential of a sum of AbstractPotentials"""
function potential(pot::Vector{<:AbstractPotential}, x::AbstractArray{T}) where {T<:Real}
    sum_pot = zeros(3)
    for i ∈ eachindex(pot)
        sum_pot .+= potential(pot[i], x)
    end
    return sum_pot
end


"""Plummer potential"""
function potential(pot::Plummer, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / sqrt(pot.b^2 + x'x)
end

"""Unitful Plummer potential"""
function potential(pot::Plummer, x::Vector{<:U.Length})
    return uconvert(u_V^2, -Gᵤ*pot.m_u / sqrt(pot.b_u^2 + x'x))
end


"""Miyamoto-Nagai disk potential"""
function potential(pot::MiyamotoNagaiDisk, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m/sqrt( x[1:2]'x[1:2] + (pot.a + sqrt(pot.b^2+x[3]^2))^2 )
end

"""Unitful MiyamotoNagaiDisk potential"""
function potential(pot::MiyamotoNagaiDisk, x::Vector{<:U.Length})
    return uconvert(u_V^2, -Gᵤ*pot.m_u/sqrt( x[1:2]'x[1:2] + (pot.a_u + sqrt(pot.b_u^2+x[3]^2))^2))
end


"""Allen and Santillan (generalized) halo"""
function potential(pot::AllenSantillanHalo, x::AbstractArray{T}) where {T<:Real}
    f(y) = 1.0 + (y/pot.a)^(pot.γ-1.0)
    r  = sqrt( x'x )
    if r < pot.Λ
        res = -G*(pot.m/pot.a)*( log(f(r)/f(pot.Λ))/(pot.γ-1.0) - (1.0-1.0/f(pot.Λ)) )
    else
        res = -G*(pot.m/r)*(pot.Λ/pot.a)^pot.γ/f(pot.Λ)
    end
    return res
end

