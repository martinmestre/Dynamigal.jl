"""Acceleration of AbstractPotentials"""


"""Acceleration"""
function acceleration(pot::AbstractPotential, x::AbstractArray{T}) where {T<:Real}
    return -1.0*gradient(y->potential(pot, y), x)[1]
end

"""Unitful acceleration"""
function acceleration(pot::UnionAbstractPotentials, x::Vector{<:U.Length})
    x = ustrip(uconvert.(u_L, x))
    return u_A*acceleration(pot, x)
end

"""Acceleration in a sum of AbstractPotentials"""
function acceleration(pot::Vector{<:AbstractPotential}, x::AbstractArray{T}) where {T<:Real}
    sum_acc = zeros(3)
    for i ∈ eachindex(pot)
        sum_acc .+= acceleration(pot[i], x)
    end
    return sum_acc
end



"""ODE for the Newtonian case"""
function ode(u,p,t)
    return SA[u[4:6]...,acceleration(p, u[1:3])...]
end