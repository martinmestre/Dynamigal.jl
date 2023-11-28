"""Acceleration of AbstractPotentials"""


"""Unitful acceleration"""
function acceleration(pot::UnionAbstractPotentials, x::Vector{<:Unitful.Length}, t::T=0𝕦.t) where {T<:Unitful.Time}
    x, t = code_units(x, t)
    return acceleration(pot, x, t)*𝕦.a
end


"""Acceleration of a sum of AbstractPotentials"""
function acceleration(pot::Vector{<:AbstractPotential}, x::AbstractArray{L}, t::T=0.0) where {L<:Real, T<:Real}
    sum_acc = zeros(3)
    for i ∈ eachindex(pot)
        sum_acc .+= acceleration(pot[i], x, t)
    end
    return sum_acc
end


"""Acceleration of single potential"""
function acceleration(pot::AbstractPotential, x::AbstractArray{L}, t::T=0.0) where {L<:Real, T<:Real}
    return -1.0*gradient(y->potential(pot, y, t), x)[1]
end


"""Complement function to be used in acceleration computation"""
selec(i::I) where {I<:Integer} = 1+3(i-1)

function complement(p::Vector{<:AbstractPotential}, x::AbstractArray{L}, i::I) where {L<:Real, I<:Integer}
    p_c = copy(p)
    x_cc = copy(x)
    deleteat!(p_c, i)
    deleteat!(x_cc, selec(i):selec(i)+2)
    x_c = Vector{Vector{Float64}}(undef,length(p_c))
    for j ∈ eachindex(x_c)
        x_c[j] = x_cc[selec(j):selec(j)+2]
    end
    return p_c, x_c
end

"""Acceleration of a system of macro particles"""
function acceleration(mps::Vector{<:AbstractMacroParticle}, x::AbstractArray{L}, t::T=0.0) where {L<:Real, T<:Real}
    p = [mps[i].pot for i ∈ eachindex(mps)] # we only use the potentials from each MacroParticle
    acc = Vector{Vector{Float64}}(undef,length(p))
    for i ∈ eachindex(p)
        y = x[selec(i):selec(i)+2]
        p_c, x_c = complement(p, x, i)
        acc[i] = sum( [acceleration(p_c, y-x_c[j]) for j ∈ eachindex(p_c)] )
        @show i acc[i]
        @show p p_c
        @show y x
        @show eachindex(p_c) eachindex(y)
    end
    return vcat(acc...)
end


