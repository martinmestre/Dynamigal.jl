"""Acceleration of AbstractPotentials"""


"""Unitful acceleration"""
function acceleration(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}; kargs...)
    x = ustrip(uconvert.(lu.l, x))
    return acceleration(pot, x; kargs...)*lu.a
end


"""Acceleration of a sum of AbstractPotentials"""
function acceleration(pot::Vector{<:AbstractPotential}, x::AbstractArray{T}) where {T<:Real}
    sum_acc = zeros(3)
    for i ∈ eachindex(pot)
        sum_acc .+= acceleration(pot[i], x; kwargs...)
    end
    return sum_acc
end


"""Acceleration of single potential"""
function acceleration(pot::AbstractPotential, x::AbstractArray{T}; kargs...) where {T<:Real}
    return -1.0*gradient(y->potential(pot, y; kwargs...), x)[1]
end







"""ODE for the Newtonian case"""
function ode(u,p,t)
    return SA[u[4:6]...,acceleration(p, u[1:3])...]
end