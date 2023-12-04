"""Potential types"""


@with_kw struct TimeDependent{T<:Real} <: AbstractPotential
    m::T
    @assert m>0 "m must be possitive"
end
TimeDependent(m::T) where {T<:Unitful.Mass} = TimeDependent( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct Kepler{T<:Real} <: AbstractPotential
    m::T
    @assert m>0 "m must be possitive"
end
Kepler(m::T) where {T<:Unitful.Mass} = Kepler( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct Plummer{T<:Real,D<:Real} <: AbstractPotential
    m::T
    a::D
    @assert m>0 & a>0 "all fields should be possitive"
end
Plummer(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct Hernquist{T<:Real,D<:Real} <: AbstractPotential
    m::T
    a::D
    @assert m>0 & a>0 "all fields should be possitive"
end
Hernquist(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Hernquist( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct MiyamotoNagaiDisk{T<:Real,D<:Real,F<:Real} <: AbstractDiskPotential
    m::T
    a::D
    b::F
    @assert m>0 & a>0 & b>0 "all fields should be possitive"
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)), ustrip(uconvert(𝕦.l, b)) )


@with_kw struct AllenSantillanHalo{T<:Real,D<:Real,F<:Real,G<:Real} <: AbstractHaloPotential
    m::T
    a::D
    Λ::F
    γ::G
    @assert m>0 & a>0 & Λ>0  "fields m, a, Λ should be possitive"
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)),
                        ustrip(uconvert(𝕦.l, Λ)),  γ )


@with_kw struct NFW{T<:Real, D<:Real, F<:Real} <: AbstractHaloPotential
    m::T  # virial mass (200)
    rₛ::D  # r_200/c
    c::F  # concentratiom
    @assert m>0 & rₛ>0 & c>0 "all fields should be possitive"
    function NFW(m::T, rₛ::D) where {T<:Real,D<:Real}
        r = (m/(200*𝕔.ρ_c*4.0/3.0*π))^(1.0/3.0)  # virial radius
        c = r/rₛ
        return New(m, rₛ, c)
    end
end
NFW(m::T, rₛ::L, c::D) where {T<:Unitful.Mass, L<:Unitful.Length, D<:Real} =
    NFW( ustrip(uconvert(𝕦.m, m)), ustrip(uconvert(𝕦.l, rₛ)), c)
