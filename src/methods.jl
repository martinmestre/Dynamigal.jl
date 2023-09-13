"""Methods"""

"""Unitless Plummer potential"""
function potential(pot::Plummer, x::AbstractArray{T}) where {T<:Real}
    return -G*pot.m / sqrt(pot.b^2 + x'x)
end

"""Unitful Plummer potential"""
function potential(pot::Plummer, x::Vector{<:U.Length})
    return uconvert(u_Pot, -Gᵤ*pot.m_u / sqrt(pot.b_u^2 + x'x))
end

"""Unitless Plummer acceleration"""
function acceleration(pot::AbstractPotential, x::AbstractArray{T}) where {T<:Real}
    return -1.0*gradient(y->potential(pot, y), x)[1]
end

"""Unitful Plummer acceleration"""
function acceleration(pot::AbstractPotential, x::Vector{<:U.Length})
    x = ustrip(uconvert.(u_L, x))
    return u_A*acceleration(pot, x)
end


function ode(u,p,t)
    return SA[u[4:6]...,acceleration(p, u[1:3])...]
end

function evolve(pot::AbstractPotential, x::Vector{<:U.Length}, v::Vector{<:U.Velocity},
                t_span::Tuple{<:U.Time, <:U.Time})
    dimₓ, dimᵥ = unit(x[begin]), unit(v[begin])
    x = ustrip(uconvert.(u_L, x))
    v = ustrip(uconvert.(u_V, v))
    t_span = uconvert.(unit(u_T),t_span)./u_T
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
    sol=solve(prob, Vern9(); abstol=5.e-8, reltol=5.e-8)
    orb = Orbit(sol.t*u_T, sol.u[1:3,:]*u_L, sol.u[4:6,:]*u_V)
    return orb
end

function test()
    plum = Plummer(10.0^11*u"Msun",10.0u"kpc")
    x₀ = [10.0, 0.0, 0.0]u"kpc"
    v₀ = [0.0,50.0,0.0]u"km/s"
    t_range = (0.0,10.0).*u_T
    @show t_range
    sol = evolve(plum, x₀, v₀, t_range)
    return sol
end