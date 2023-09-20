"""Evolution functions"""

"""Evolution of a an initial condition in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, x::Vector{<:U.Length}, v::Vector{<:U.Velocity},
    t_span::Tuple{<:U.Time, <:U.Time})
    dimₓ, dimᵥ = unit(x[begin]), unit(v[begin])
    x = ustrip(uconvert.(u_L, x))
    v = ustrip(uconvert.(u_V, v))
    t_span = uconvert.(unit(u_T),t_span)./u_T
    p = pot
    u₀ = SA[x...,v...]
    prob = ODEProblem(ode, u₀, t_span, p)
    sol=solve(prob, Vern9(); abstol=5.e-8, reltol=5.e-8)
    orb = Orbit(sol.t*u_T, sol[1:3,:]*u_L, sol[4:6,:]*u_V)
    return orb
end

"""Evolution of a TestParticle in an AbstractPotential"""
function evolve(pot::UnionAbstractPotentials, p::AbstractParticle, t_span::Tuple{<:U.Time, <:U.Time})
    x = p.x_u
    v = p.v_u
    return evolve(pot, x, v, t_span)
end



