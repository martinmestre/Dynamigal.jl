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


@with_kw struct Hernquist{T<:Real,D<:Real} <: AbstractSphericalStaticPotential # mass is not defined yet but acceleration is defined so it dispatchs correctly.
    m::T
    a::D
    @assert m>0 && a>0 "all fields should be possitive"
end
Hernquist(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Hernquist( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct Kepler{T<:Real} <: AbstractStaticPotential # do not assign spherical as mass(r) is not defined yet.
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
    ρ₀::D  # central density = (1/3) * m / volume(a)
    𝔸::D  # NFW constant: 𝔸 = A_nfw = f_nfw(c)
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

function NFW_from_m_c(m::T, c::F, cosmos::C=𝕔) where {T<:Real,F<:Real,C<:AbstractConfig}
    @assert m>0 && c>0  "all fields should be possitive"
    𝔸 = f_nfw(c)
    m_v = 𝔸*m
    r_v = r_vir_nfw(m_v, cosmos)
    a = r_v/c
    ρ₀ = m / (4π*a^3)
    return NFW(m, a, c, m_v, r_v, ρ₀, 𝔸, cosmos)
end
NFW_from_m_c(m::M, c::L, cosmos::C=𝕔) where {M<:Unitful.Mass, L<:Real, C<:AbstractConfig} =
     NFW( ustrip(uconvert(𝕦.m, m)), c, cosmos)
NFW_from_m_c(; m::M, c::L, cosmos::C=𝕔) where {M,L,C} = NFW_from_m_c(m, c, cosmos)

function NFW_from_mv_a(m_v::T, a::F, cosmos::C=𝕔) where {T<:Real,F<:Real,C<:AbstractConfig}
    @assert m_v>0 && a>0  "all fields should be possitive"
    r_v = r_vir_nfw(m_v, cosmos)  # virial radius
    c = r_v/a
    𝔸 = f_nfw(c)
    m = m_v/𝔸
    ρ₀ = m / (4π*a^3)
    return NFW(m, a, c, m_v, r_v, ρ₀, 𝔸, cosmos)
end
NFW_from_mv_a(m_v::M, a::L, cosmos::C=𝕔) where {M<:Unitful.Mass, L<:Unitful.Length, C<:AbstractConfig} =
     NFW_from_mv_a( ustrip(uconvert(𝕦.m, m_v)), ustrip(uconvert(𝕦.l, a)), cosmos)
NFW_from_mv_a(; m_v::M, a::L, cosmos::C=𝕔) where {M,L,C} = NFW_from_mv_a(m_v, c, cosmos)

function NFW_from_mv_c(m_v::T, c::F, cosmos::C=𝕔) where {T<:Real,F<:Real,C<:AbstractConfig}
    @assert m_v>0 && c>0  "all fields should be possitive"
    r_v = r_vir_nfw(m_v, cosmos)  # virial radius
    a = r_v/c
    𝔸 = f_nfw(c)
    m = m_v/𝔸
    ρ₀ = m / (4π*a^3)
    return NFW(m, a, c, m_v, r_v, ρ₀, 𝔸, cosmos)
end
NFW_from_mv_c(m_v::M, c::L, cosmos::C=𝕔) where {M<:Unitful.Mass, L<:Real, C<:AbstractConfig} =
     NFW( ustrip(uconvert(𝕦.m, m_v)),  c, cosmos)
NFW_from_mv_c(; m_v::M, c::L, cosmos::C=𝕔) where {M,L,C} = NFW_from_mv_c(m_v, c, cosmos)

"""General NFW with symbol arguments."""
function NFW(s₁::Symbol, s₂::Symbol, q₁::T, q₂::F, cosmos::C=𝕔) where {T, F, C<:AbstractConfig}
    if s₁==:m && s₂==:a
        return NFW(q₁, q₂, cosmos)
    elseif s₁==:m && s₂==:c
        return NFW_from_m_c(q₁, q₂, cosmos)
    elseif s₁==:m_v && s₂==:a
        return NFW_from_mv_a(q₁, q₂, cosmos)
    elseif s₁==:m_v && s₂==:c
        return NFW_from_mv_c(q₁, q₂, cosmos)
    else
        error("Symbols should be any of the following combinations
        :m, :a,
        :m, :c
        :m_v, :a
        :m_v, :c
        ")
    end
end


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


struct PowerLawCutoff{T<:Real,D<:Real,R<:Real,G<:Real} <: AbstractSphericalStaticPotential
    m::T # total mass
    α::D # power-law index
    c::R # cutoff radius
    β::G # auxiliary constant
    Γ_β::G # Γ(β)
    𝔸::G  # aux constant
end
function PowerLawCutoff(m::T,α::D,c::R) where {T<:Real,D<:Real,R<:Real}
    @assert m>0 && α>=0 && α<3 && c>0 "Restrictions: m > 0, 0 ≤ α < 3, c > 0"
    β = 0.5*(3-α) # auxiliary constant
    Γ_β = gamma(β)
    𝔸 = (m/2π)*c^(α-3)/Γ_β
    return PowerLawCutoff(m, α, c, β, Γ_β, 𝔸)
end
PowerLawCutoff(m::T, α::D, c::R) where {T<:Unitful.Mass, D<:Real, R<:Unitful.Length} =
     PowerLawCutoff( ustrip(uconvert(𝕦.m, m)), α, ustrip(uconvert(𝕦.l, c)) )
function PowerLawCutoff(; m::T, α::D, c::R) where {T,D,R}
    if typeof(m)<:Unitful.Mass && typeof(c)<:Unitful.Length
        return PowerLawCutoff(m, α, c)
    elseif typeof(m)<:Real && typeof(c)<:Real
        return PowerLawCutoff(m, α, c)
    else
        error("m and c shoud be both quantities or both numbers")
    end
end

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


