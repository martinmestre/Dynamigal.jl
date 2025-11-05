"""Potential types"""

@with_kw struct AllenSantillanHalo{T<:Real,D<:Real,F<:Real,G<:Real} <: AbstractSphericalStaticPotential
    m::T
    a::D
    Λ::F
    γ::G
    @assert m>0 && a>0 && Λ>0  "fields m, a, Λ should be possitive"
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)),
                        ustrip(uconvert(𝕦.l, Λ)),  γ )


@with_kw struct Hernquist{T<:Real,D<:Real} <: AbstractSphericalStaticPotential
    m::T
    a::D
    @assert m>0 && a>0 "all fields should be possitive"
end
Hernquist(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Hernquist( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct Kepler{T<:Real} <: AbstractSphericalStaticPotential
    m::T
    @assert m>0 "m must be possitive"
end
Kepler(m::T) where {T<:Unitful.Mass} = Kepler( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct MiyamotoNagaiDisk{T<:Real,D<:Real,F<:Real} <: AbstractStaticPotential
    m::T
    a::D
    b::F
    @assert m>0 && a>0 && b>0 "all fields should be possitive"
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)), ustrip(uconvert(𝕦.l, b)) )


# NFW
f_nfw(x::T) where {T<:Real} = log(1+x)-x/(1+x)

function r_vir_nfw(m_v::T, cosmos::C=𝕔)  where {T<:Real,C<:AbstractConfig}
    ρ_c = ustrip( uconvert(𝕦.m/𝕦.l^3, cosmos.ρ_c))
    r_v = (m_v/(200*ρ_c*4.0/3.0*π))^(1.0/3.0)
    return r_v
end
r_vir_nfw(m_v::M, cosmos::C=𝕔) where {M<:Unitful.Mass,C<:AbstractConfig} =
    r_vir_nfw(adimensional(m_v), cosmos)

function concentration(m::T, a::F, cosmos::C) where {T<:Real,F<:Real,C<:AbstractConfig}
    ρ_c = ustrip( uconvert(𝕦.m/𝕦.l^3, cosmos.ρ_c))
    𝔹 = (4π/3)*200*ρ_c*a^3/m
    g(x) = f_nfw(x)/x^3 - 𝔹
    D(f)= x->gradient(y->f(y),x)[1]
    return find_zero((g,D(g)),  [1.0e-6,100.0], Roots.Brent())
end

"""
    NFW struct
"""
struct NFW{T<:Real, F<:Real, D<:Real, C<:AbstractConfig} <: AbstractSphericalStaticPotential
    m::T  # scale mass
    a::F  # scale radius
    c::D  # concentration: c=r_v/a
    m_v::D # virial mass
    r_v::D  # virial radius
    ρ₀::D  # central density
    𝔸::D  # NFW constant: A_nfw = f_nfw(c)
    cosmos::C  # configuration struct
end

function NFW(m::T, a::F, cosmos::C=𝕔) where {T<:Real,F<:Real,C<:AbstractConfig}
    @assert m>0 && a>0  "all fields should be possitive"
    c = concentration(m, a, cosmos)
    𝔸 = f_nfw(c)
    m_v = 𝔸*m
    r_v = c*a
    ρ₀ = m / (4π*a^3)
    return NFW(m, a, c, m_v, r_v, ρ₀, 𝔸, cosmos)
end
NFW(m::M, a::L, cosmos::C=𝕔) where {M<:Unitful.Mass, L<:Unitful.Length, C<:AbstractConfig} =
     NFW( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)), cosmos)
NFW(; m::M, a::L, cosmos::C=𝕔) where {M,L,C} = NFW(m, a, cosmos)

"""NFW from virial mass and concentration. Only arguments, not keywargs. m_v = Unitful Q."""
function NFW(m_v::T, c::F, cosmos::C=𝕔) where {T<:Unitful.Mass, F<:Real, C<:AbstractConfig}
    m_v = adimensional(m_v)
    @assert m_v>0 && c>0  "all fields should be possitive"
    r_v = r_vir_nfw(m_v, cosmos)  # virial radius
    a = r_v/c
    𝔸 = f_nfw(c)
    m = m_v/𝔸
    ρ₀ = m / (4π*a^3)
    return NFW(m, a, c, m_v, r_v, ρ₀, 𝔸, cosmos)
end


# @with_kw struct NFW{T<:Real, F<:Real, D<:Real, C<:AbstractConfig} <: AbstractSphericalStaticPotential
#     @assert m_v>0 && a>0  "all fields should be possitive"
#     m_v::T  # virial mass: M(r_v)
#     a::F  # scale radius: a=r_v/c
#     𝕔_l::C = 𝕔
#     r_v::D = r_vir_nfw(m_v; 𝕔_l=𝕔_l) # virial radius
#     c::D = r_v/a # concentration: c=r_v/a
#     𝔸::D = f_nfw(c)
#     m::D = m_v/𝔸 # scale mass = m_v / f_nfw(c)
#     ρ₀::D = m / (4π*a^3) # central density
# end
# NFW(m_v::T, a::F; 𝕔_l=𝕔) where {T,F} = NFW(; m_v=m_v, a=a, 𝕔_l=𝕔_l)
# NFW(m_v::M, a::L) where {M<:Unitful.Mass, L<:Unitful.Length} =
#     NFW( ustrip(uconvert(𝕦.m, m_v)),  ustrip(uconvert(𝕦.l, a)))

# function NFW(; m::T, a::F, 𝕔_l=𝕔) where {T<:Real,F<:Real}
#     c = concentration(m, a, 𝕔_l)
#     𝔸 = f_nfw(c)
#     m_v = 𝔸*m
#     r_v = c*a
#     ρ₀ = m / (4π*a^3)
#     return NFW(m_v, a, 𝕔_l, r_v, c, 𝔸, m, ρ₀)
# end

# function NFW(m_v::M, c::T; 𝕔_l=𝕔) where {M<:Unitful.Mass, T<:Real}
#     m_v = adimensional(m_v)
#     @assert m_v>0 && c>0  "all fields should be possitive"
#     r_v = r_vir_nfw(m_v; 𝕔_l=𝕔_l)  # virial radius
#     a = r_v/c
#     𝔸 = f_nfw(c)
#     m = m_v/𝔸
#     ρ₀ = m / (4π*a^3)
#     return NFW(m_v, a, 𝕔_l, r_v, c, 𝔸, m, ρ₀)
# end



"""Time dependent potentials"""
@with_kw struct OscillatoryKepler{T<:Real, D<:Real} <: AbstractPotential
    m::T
    τ::D
    @assert (m>0 && τ::D) "all fields should be possitive"
end
OscillatoryKepler(m::T, τ::D) where {T<:Unitful.Mass, D<:Unitful.Time} = OscillatoryKepler( ustrip(uconvert(𝕦.m, m)), ustrip(uconvert(𝕦.t, τ) ) )
time_dependence(::Type{<:OscillatoryKepler}) = TimeDependent()


@with_kw struct Plummer{T<:Real,D<:Real} <: AbstractSphericalStaticPotential
    m::T
    a::D
    @assert (m>0 && a>0) "all fields should be possitive"
end
Plummer(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct PowerLawCutoff{T<:Real,D<:Real,R<:Real,G<:Real} <: AbstractSphericalStaticPotential
    m::T # total mass
    α::D # power-law index
    c::R # cutoff radius
    β::G = 0.5*(3-α) # auxiliary constant
    𝔸::G = (m/2π)*c^(α-3)/gamma(β)
    @assert m>0 && α>=0 && α<3 && c>0 "all fields should be possitive"
end
PowerLawCutoff(m::T, α::D, c::R) where {T<:Unitful.Mass, D<:Real, R<:Unitful.Length} =
    PowerLawCutoff( ustrip(uconvert(𝕦.m, m)), α, ustrip(uconvert(𝕦.l, c)) )


"""Composite types"""
struct CompositePotential{P <: NTuple{N, AbstractPotential} where N} <: AbstractPotential
    potentials::P
end

function CompositePotential(potentials::NTuple{N, T}) where {N, T <: AbstractPotential}
    if N == 1
        throw(ArgumentError("CompositePotential requires at least 2 elements, got $N"))
    end
    return CompositePotential{typeof(potentials)}(potentials)
end
CompositePotential(p...) = CompositePotential(p)


