
"""ODE for the Newtonian case: test particle in a potential."""
function ode(u::AbstractArray{<:Real}, p::UnionAbstractPotentials, t::T) where {T<:Real}
    return SA[u[siss]..., acceleration(p, u[sis], t)...]
end

"""ODE for the Newtonian case: system of macro particles."""
function ode(u::AbstractArray{<:Real}, p::Vector{<:AbstractMacroParticle}, t::T) where {T<:Real}
    n = Integer(length(u)/2)
    return SA[u[n+1:end]..., acceleration(p, u[begin:n], t)...]
end