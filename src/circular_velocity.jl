"""Circular velocity"""

function circular_velocity(pot::UnionAbstractPotentials, x::AbstractArray{L}, t::T=0.0) where {L<:Real, T<:Real}
    a = acceleration(pot, x, t)
    return sqrt( sqrt(x'x)*sqrt(a'a) )
end

function circular_velocity(pot::UnionAbstractPotentials, r::L, t::T=0.0) where {L<:Real, T<:Real}
    try sqrt(r)
        x = [r, 0., 0.]
        return circular_velocity(pot, x, t)
    catch err
        rethow(err)
    end
end

function circular_velocity(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, t::T=0𝕦.t) where {T<:Unitful.Time}
    x, t = code_units(x, t)
    return circular_velocity(pot, x, t)*𝕦.v
end

function circular_velocity(pot::UnionAbstractPotentials, r::L, t::T=0𝕦.t) where {L<:Unitful.Length, T<:Unitful.Time}
    x = [r, 0.0𝕦.l, 0.0𝕦.l]
    return circular_velocity(pot, x, t)
end