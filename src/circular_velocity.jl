"""Circular velocity"""

function circular_velocity(pot::P, x::AbstractArray{L}, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    a = acceleration(pot, x, t)
    return sqrt( sqrt(x'x)*sqrt(a'a) )
end

# function circular_velocity(pot::P, r::L, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
#     try sqrt(r)
#         x = [r, 0., 0.]
#         return circular_velocity(pot, x, t)
#     catch err
#         rethrow(err)
#     end
# end
function circular_velocity(pot::P, r::L, t::T=0.0) where {P<:AbstractPotential, L<:Real, T<:Real}
    x = [r, 0., 0.]
    return circular_velocity(pot, x, t)
end

function circular_velocity(pot::P, x::Vector{<:Unitful.Length}, t::T=0𝕦.t) where {P<:AbstractPotential,T<:Unitful.Time}
    x, t = code_units(x, t)
    return circular_velocity(pot, x, t)*𝕦.v
end

function circular_velocity(pot::P, r::L, t::T=0𝕦.t) where {P<:AbstractPotential, L<:Unitful.Length, T<:Unitful.Time}
    x = [r, 0.0𝕦.l, 0.0𝕦.l]
    return circular_velocity(pot, x, t)
end