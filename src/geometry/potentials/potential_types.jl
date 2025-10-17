"""Potential types"""


@with_kw struct Kepler{T<:Real} <: AbstractStaticPotential
    m::T
    @assert m>0 "m must be possitive"
end
Kepler(m::T) where {T<:Unitful.Mass} = Kepler( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct Plummer{T<:Real,D<:Real} <: AbstractStaticPotential
    m::T
    a::D
    @assert (m>0 && a>0) "all fields should be possitive"
end
Plummer(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct Hernquist{T<:Real,D<:Real} <: AbstractStaticPotential
    m::T
    a::D
    @assert m>0 && a>0 "all fields should be possitive"
end
Hernquist(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Hernquist( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct PowerLawCutoff{T<:Real,D<:Real,R<:Real} <: AbstractStaticPotential
    m::T # total mass
    α::D # power-law index
    c::R # cutoff radius
    @assert m>0 && a>0 && c>R "all fields should be possitive"
end
PowerLawCutoff(m::T, α::D, c::R) where {T<:Unitful.Mass, D<:Real, T<:Unitful.Length} =
    PowerLawCutoff( ustrip(uconvert(𝕦.m, m)), α, ustrip(uconvert(𝕦.l, a)) )


@with_kw struct MiyamotoNagaiDisk{T<:Real,D<:Real,F<:Real} <: AbstractStaticPotential
    m::T
    a::D
    b::F
    @assert m>0 && a>0 && b>0 "all fields should be possitive"
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)), ustrip(uconvert(𝕦.l, b)) )


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



# NFW
f_nfw(x::T) where {T<:Real} = log(1+x)-x/(1+x)

function r_vir_nfw(m; 𝕔=𝕔)
    ρ = ustrip( uconvert(𝕦.m/𝕦.l^3, 𝕔.ρ_c))
    r = (m/(200*ρ*4.0/3.0*π))^(1.0/3.0)
    return r
end
r_vir_nfw(m::M; 𝕔=𝕔) where {M<:Unitful.Mass} = r_vir_nfw(adimensional(m); 𝕔=𝕔)


@with_kw struct NFW{T<:Real, F<:Real, D<:Real, C<:AbstractConfig} <: AbstractStaticPotential
    @assert m>0 && a>0  "all fields should be possitive"
    m::T  # virial mass: M(r)
    a::F  # scale radius: a=r/c
    𝕔::C = 𝕔
    r::D = r_vir_nfw(m; 𝕔=𝕔) # virial radius
    c::D = r/a # concentration: c=r/a
    𝔸::D = f_nfw(c)
end
NFW(m::T, a::F; 𝕔=𝕔) where {T,F} = NFW(; m=m, a=a, 𝕔=𝕔)
NFW(m::M, a::L) where {M<:Unitful.Mass, L<:Unitful.Length} =
    NFW( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)))


function NFW(m::M, c::T; 𝕔=𝕔) where {M<:Unitful.Mass, T<:Real}
    m = adimensional(m)
    @assert m>0 && c>0  "all fields should be possitive"
    r = r_vir_nfw(m; 𝕔=𝕔)  # virial radius
    a = r/c
    𝔸 = f_nfw(c)
    return NFW(m, a, 𝕔, r, c, 𝔸)
end

function concentration(p::NFW; 𝕔=𝕔)
    ρ = ustrip( uconvert(𝕦.m/𝕦.l^3, 𝕔.ρ_c))
    r = (p.m/(200*ρ*4.0/3.0*π))^(1.0/3.0)  # virial radius
    return r/p.a
end


"""Time dependent potentials"""
@with_kw struct OscillatoryKepler{T<:Real, D<:Real} <: AbstractPotential
    m::T
    τ::D
    @assert (m>0 && τ::D) "all fields should be possitive"
end
OscillatoryKepler(m::T, τ::D) where {T<:Unitful.Mass, D<:Unitful.Time} = OscillatoryKepler( ustrip(uconvert(𝕦.m, m)), ustrip(uconvert(𝕦.t, τ) ) )
time_dependence(::Type{<:OscillatoryKepler}) = TimeDependent()



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


"""
Customized potential type constructors
"""

"""MilkyWayBovy2014
MWPotential2014 = [
    PowerSphericalPotentialwCutoff(normalize=0.05, alpha=1.8, rc=1.9 / 8.0),
    MiyamotoNagaiPotential(a=3.0 / 8.0, b=0.28 / 8.0, normalize=0.6),
    NFWPotentia
    l(a=2.0, normalize=0.35),
]
https://github.com/jobovy/galpy/blob/b7c3bf055880d21f4b250981acfdd5f0b4f5db09/galpy/potential/mwpotentials.py#L24
"""
function MilkyWayBovy2014()
    bulge = PowerLawCutoff(m=, a=)
    disk = MiyamotoNagaiDisk(m= , a=3.0, b=0.28)
    halo = NFW(m= , a=)
    return CompositePotential(bulge,disk,halo)
end