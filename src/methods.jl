"""Methods"""

"""Adimensional Plummer potential"""
function potential(pot::Plummer, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / sqrt(pot.b^2 + x'x)
end

"""Unitful Plummer potential"""
function potential(pot::Plummer, x::Vector{<:U.Length})
    return uconvert(u_Pot, -Gᵤ*pot.m_u / sqrt(pot.b_u^2 + x'x))
end

"""Adimensional Plummer acceleration"""
function acceleration(pot::AbstractPotential, x::AbstractArray{T}) where {T<:Real}
    return -1.0*gradient(y->potential(pot, y), x)[1]
end

"""Unitful Plummer acceleration"""
function acceleration(pot::AbstractPotential, x::Vector{<:U.Length})
    x = [uconvert(u_L, x[i]).val for i ∈ eachindex(x)]
    return u_A*acceleration(pot, x)
end


function ode(u,p,t)
    return SA[u[4:6]...,acceleration(p[1], u[1:3])...]
end

function test()
    plum = Plummer(10.0^11*u"Msun",10.0u"kpc")
    p = [plum]
    x₀ = [10.0, 0.0, 0.0]
    v₀ = [0.0,50.0,0.0]
    u₀ = SA[x₀...,v₀...]
    tspan = (0.0, 10.0)
    prob = ODEProblem(ode, u₀, tspan, p)
    sol=solve(prob, Tsit5(); abstol=5.e-8, reltol=5.e-8)
    return sol
end