
"""ODE for the Newtonian case: test particle in a potential."""
function ode(u::AbstractArray{<:Real}, p::P, t::T) where {P<:AbstractPotential, T<:Real}
    return SA[u[siss]..., acceleration(p, u[sis], t)...]
end

"""ODE for the Newtonian case: system of macro particles."""
function ode(u::AbstractArray{<:Real}, p::Vector{<:AbstractMacroParticle}, t::T) where {T<:Real}
    n = Integer(length(u)/2)
    return SA[u[n+1:end]..., acceleration(p, u[begin:n], t)...]
end


"""Sugerencia ..."""
function ode!(du, u, p::Vector{<:AbstractMacroParticle}, t::T) where {T<:Real}
    n = length(p)
    @inbounds for i in 1:n
        # posiciones y velocidades
        xi = @view u[selec(i):selec(i)+2]
        vi = @view u[n*3+selec(i):n*3+selec(i)+2]

        # acumulador de aceleración
        acc_i = MVector{3,Float64}(0.0,0.0,0.0)

        @inbounds for j in 1:n
            if j != i
                xj = @view u[selec(j):selec(j)+2]
                acc_i .+= acceleration(p[j].pot, xi - xj, t)
            end
        end

        # escribir resultado en du
        du[selec(i):selec(i)+2] .= vi          # derivada de posición = velocidad
        du[n*3+selec(i):n*3+selec(i)+2] .= acc_i  # derivada de velocidad = aceleración
    end
    return nothing
end
