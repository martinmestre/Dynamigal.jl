"""Potential types"""


@with_kw struct TimeDependent{T<:Real} <: AbstractPotential
    m::T
end
TimeDependent(m::T) where {T<:Unitful.Mass} = TimeDependent( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct PointMass{T<:Real} <: AbstractPotential
    m::T
end
PointMass(m::T) where {T<:Unitful.Mass} = PointMass( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct Plummer{T<:Real,D<:Real} <: AbstractPotential
    m::T
    b::D
end
Plummer(m::T, b::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, b)) )


@with_kw struct MiyamotoNagaiDisk{T<:Real,D<:Real,F<:Real} <: AbstractDiskPotential
    m::T
    a::D
    b::F
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)), ustrip(uconvert(𝕦.l, b)) )


@with_kw struct AllenSantillanHalo{T<:Real,D<:Real,F<:Real,G<:Real} <: AbstractDiskPotential
    m::T
    a::D
    Λ::F
    γ::G
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)),
                        ustrip(uconvert(𝕦.l, Λ)),  γ )
