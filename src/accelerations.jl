"""Acceleration of AbstractPotentials"""


"""Unitful acceleration"""
function acceleration(pot::P, x::Vector{<:Unitful.Length}, t::T) where {P<:AbstractPotential, T<:Unitful.Time}
    x, t = adimensional(x, t)
    return acceleration(pot, x, t)*𝕦.a
end
function acceleration(pot::P, x::Vector{<:Unitful.Length}) where {P<:AbstractPotential}
    x = adimensional(x)
    return acceleration(pot, x)*𝕦.a
end


"""Acceleration of a sum of AbstractPotentials"""
function acceleration(pot::CompositePotential, x::AbstractArray{L}, t::T=0.0) where {L<:Real, T<:Real}
    sum_acc = zeros(L, 3)
    for p ∈ pot
        sum_acc .+= acceleration(p, x, t)
    end
    return sum_acc
end

"""
    acceleration(pot::P, x::AbstractArray{L}, t::T) where {P<:AbstractStaticPotential, L<:Real, T<:Real}
Bridge function for static potentials
"""
function acceleration(pot::P, x::AbstractArray{L}, t::T) where {P<:AbstractStaticPotential, L<:Real, T<:Real}
    return acceleration(pot, x)
end

"""Acceleration of single potential. Static and time-dependent method."""
function acceleration(pot::P, x::AbstractArray{L}, t::T) where {P<:AbstractPotential, L<:Real, T<:Real}
    return -gradient(y->potential(pot, y, t), x)[1]
end
function acceleration(pot::P, x::AbstractArray{L}) where {P<:AbstractPotential, L<:Real}
    return -gradient(y->potential(pot, y), x)[1]
end

"""Complement function to be used in acceleration computation"""
selec(i::I) where {I<:Integer} = 1+3(i-1)

function complement(p::Vector{<:AbstractPotential}, x::AbstractArray{L}, i::I) where {L<:Real, I<:Integer}
    ρ = copy(p)
    x_c = copy(x)
    deleteat!(ρ, i)
    deleteat!(x_c, selec(i):selec(i)+2)
    χ = Vector{Vector{Float64}}(undef,length(ρ))
    for j ∈ eachindex(χ)
        χ[j] = x_c[selec(j):selec(j)+2]
    end
    return ρ, χ
end

# function complement(p::Vector{<:AbstractPotential}, x::AbstractArray{L}, i::I) where {L<:Real, I<:Integer}
#     n = length(p)
#     ρ = [p[j] for j in 1:n if j != i]
#     χ = [SVector{3,L}(
#             x[selec(j)],
#             x[selec(j)+1],
#             x[selec(j)+2]
#         ) for j in 1:n if j != i]

#     return ρ, χ
# end

"""Acceleration of a system of macro particles"""
function acceleration(mps::Vector{<:AbstractMacroParticle}, x::AbstractArray{L}, t::T=0.0) where {L<:Real, T<:Real}
    p = [mps[i].pot for i ∈ eachindex(mps)] # we only use the potentials from each MacroParticle
    acc = Vector{Vector{Float64}}(undef,length(p))
    for i ∈ eachindex(p)
        y = x[selec(i):selec(i)+2]
        ρ, χ = complement(p, x, i)
        acc[i] = sum( [acceleration(ρ, y-χ[j]) for j ∈ eachindex(ρ)] )
    end
    return vcat(acc...)
end


"""Analytical accelerations"""

"""NFW halo acceleration"""
function acceleration(pot::NFW, x::AbstractArray{T}, t::T=0.0) where {T<:Real}
    @unpack m, a, 𝔸 = pot
    r = sqrt(x'x)
    𝕗 = -G*m/𝔸*f_nfw(r/a)/r^2
    return 𝕗*x/r
end