"""Potential types"""


@with_kw struct TimeDependent{T} <: AbstractPotential
    m::T
end
TimeDependent(m::T) where {T<:Unitful.Mass} = TimeDependent( ustrip(uconvert(lu.m, m)) )


@with_kw struct PointMass{T} <: AbstractPotential
    m::T
end
PointMass(m::T) where {T<:Unitful.Mass} = PointMass( ustrip(uconvert(lu.m, m)) )


@with_kw struct Plummer{T,D} <: AbstractPotential
    m::T
    b::D
end
Plummer(m::T, b::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(lu.m, m)),  ustrip(uconvert(lu.l, b)) )


@with_kw struct MiyamotoNagaiDisk{T,D,F} <: AbstractDiskPotential
    m::T
    a::D
    b::F
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(lu.m, m)),  ustrip(uconvert(lu.l, a)), ustrip(uconvert(lu.l, b)) )


@with_kw struct AllenSantillanHalo{T,D,F,G} <: AbstractDiskPotential
    m::T
    a::D
    Λ::F
    γ::G
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(lu.m, m)),  ustrip(uconvert(lu.l, a)),
                        ustrip(uconvert(lu.l, Λ)),  γ )
