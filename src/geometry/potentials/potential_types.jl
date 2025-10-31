"""Potential types"""

@with_kw struct AllenSantillanHalo{T<:Real,D<:Real,F<:Real,G<:Real} <: AbstractStaticPotential
    m::T
    a::D
    Λ::F
    γ::G
    @assert m>0 && a>0 && Λ>0  "fields m, a, Λ should be possitive"
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)),
                        ustrip(uconvert(𝕦.l, Λ)),  γ )


@with_kw struct Hernquist{T<:Real,D<:Real} <: AbstractStaticPotential
    m::T
    a::D
    @assert m>0 && a>0 "all fields should be possitive"
end
Hernquist(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Hernquist( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct Kepler{T<:Real} <: AbstractStaticPotential
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

function r_vir_nfw(m; 𝕔_l=𝕔)
    ρ_c = ustrip( uconvert(𝕦.m/𝕦.l^3, 𝕔_l.ρ_c))
    r_v = (m/(200*ρ_c*4.0/3.0*π))^(1.0/3.0)
    return r_v
end
r_vir_nfw(m::M; 𝕔_l=𝕔) where {M<:Unitful.Mass} = r_vir_nfw(adimensional(m); 𝕔_l=𝕔_l)

function concentration(pot::NFW)
    return r_vir_nfw(pot.m; 𝕔_l=pot.𝕔_l)/pot.a
end
"""
    NFW struct
to do: define method to enter m = scale mass and a=length scale, needs root finding for concentration "c".
"""
@with_kw struct NFW{T<:Real, F<:Real, D<:Real, C<:AbstractConfig} <: AbstractStaticPotential
    @assert m_v>0 && a>0  "all fields should be possitive"
    m_v::T  # virial mass: M(r_v)
    a::F  # scale radius: a=r_v/c
    𝕔_l::C = 𝕔
    r_v::D = r_vir_nfw(m_v; 𝕔_l=𝕔_l) # virial radius
    c::D = r_v/a # concentration: c=r_v/a
    𝔸::D = f_nfw(c)
    m::D = m_v/𝔸 # scale mass = m_v / f_nfw(c)
    ρ₀::D = m / (4π*a^3) # central density
end
NFW(m_v::T, a::F; 𝕔_l=𝕔) where {T,F} = NFW(; m_v=m_v, a=a, 𝕔_l=𝕔_l)
NFW(m_v::M, a::L) where {M<:Unitful.Mass, L<:Unitful.Length} =
    NFW( ustrip(uconvert(𝕦.m, m_v)),  ustrip(uconvert(𝕦.l, a)))


function NFW(m_v::M, c::T; 𝕔_l=𝕔) where {M<:Unitful.Mass, T<:Real}
    m_v = adimensional(m_v)
    @assert m_v>0 && c>0  "all fields should be possitive"
    r_v = r_vir_nfw(m_v; 𝕔_l=𝕔_l)  # virial radius
    a = r_v/c
    𝔸 = f_nfw(c)
    m = m_v/𝔸
    ρ₀ = m / (4π*a^3)
    return NFW(m_v, a, 𝕔_l, r_v, c, 𝔸, m, ρ₀)
end


"""Time dependent potentials"""
@with_kw struct OscillatoryKepler{T<:Real, D<:Real} <: AbstractPotential
    m::T
    τ::D
    @assert (m>0 && τ::D) "all fields should be possitive"
end
OscillatoryKepler(m::T, τ::D) where {T<:Unitful.Mass, D<:Unitful.Time} = OscillatoryKepler( ustrip(uconvert(𝕦.m, m)), ustrip(uconvert(𝕦.t, τ) ) )
time_dependence(::Type{<:OscillatoryKepler}) = TimeDependent()


@with_kw struct Plummer{T<:Real,D<:Real} <: AbstractStaticPotential
    m::T
    a::D
    @assert (m>0 && a>0) "all fields should be possitive"
end
Plummer(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct PowerLawCutoff{T<:Real,D<:Real,R<:Real} <: AbstractStaticPotential
    m::T # total mass
    α::D # power-law index
    c::R # cutoff radius
    @assert m>0 && α>=0 && α<3 && c>0 "all fields should be possitive"
end
PowerLawCutoff(m::T, α::D, c::R) where {T<:Unitful.Mass, D<:Real, R<:Unitful.Length} =
    PowerLawCutoff( ustrip(uconvert(𝕦.m, m)), α, ustrip(uconvert(𝕦.l, a)) )


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


