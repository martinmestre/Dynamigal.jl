"""Potential functions"""

"""Plummer potential"""
function potential(pot::Plummer, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / sqrt(pot.b^2 + x'x)
end

"""Unitful Plummer potential"""
function potential(pot::Plummer, x::Vector{<:U.Length})
    return uconvert(u_Pot, -Gᵤ*pot.m_u / sqrt(pot.b_u^2 + x'x))
end


"""Miyamoto-Nagai disk potential"""
function potential(pot::MiyamotoNagaiDisk, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m/sqrt( x[1:2]'x[1:2] + (pot.a + sqrt(pot.b^2+x[3]^2))^2 )
end

"""Unitful MiyamotoNagaiDisk potential"""
function potential(pot::MiyamotoNagaiDisk, x::Vector{<:U.Length})
    return uconvert(u_Pot, -Gᵤ*pot.m_u/sqrt( x[1:2]'x[1:2] + (pot.a_u + sqrt(pot.b_u^2+x[3]^2))^2))
end