"""Acceleration of AbstractPotentials"""


"""Unitful acceleration"""
function acceleration(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, t::T=nothing) where {T<:Union{Unitful.Time,Nothing}}
    x, t = code_units(x, t)
    return acceleration(pot, x, t)*𝕦.a
end


"""Acceleration of a sum of AbstractPotentials"""
function acceleration(pot::Vector{<:AbstractPotential}, x::AbstractArray{L}, t::T=nothing) where {L<:Real, T<:Union{Real,Nothing}}
    sum_acc = zeros(3)
    for i ∈ eachindex(pot)
        sum_acc .+= acceleration(pot[i], x, t)
    end
    return sum_acc
end


"""Acceleration of single potential"""
function acceleration(pot::AbstractPotential, x::AbstractArray{L}, t::T=nothing) where {L<:Real, T<:Union{Real,Nothing}}
    return -1.0*gradient(y->potential(pot, y, t), x)[1]
end




"""ODE for the Newtonian case"""
function ode(u,p,t)
    return SA[u[4:6]..., acceleration(p, u[1:3], t)...]
end