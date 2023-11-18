"""Potential types"""

@with_kw struct PointMass{T} <: AbstractPotential
    m::T
end
    PointMass(m::T) where {T<:U.Mass} = PointMass( ustrip(uconvert(u_M, m)) )


@with_kw struct Plummer{T,D} <: AbstractPotential
    m::T
    b::D
end
Plummer(m::T, b::D) where {T<:U.Mass, D<:U.Length} =
    Plummer( ustrip(uconvert(u_M, m)),  ustrip(uconvert(u_L, b)) )


@with_kw struct MiyamotoNagaiDisk{T,D,F} <: AbstractDiskPotential
    m::T
    a::D
    b::F
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:U.Mass, D<:U.Length, F<:U.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(u_M, m)),  ustrip(uconvert(u_L, a)), ustrip(uconvert(u_L, b)) )


@with_kw struct AllenSantillanHalo{T,D,F,G} <: AbstractDiskPotentia
    m::T
    a::D
    Λ::F
    γ::G
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:U.Mass, D<:U.Length, F<:U.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(u_M, m)),  ustrip(uconvert(u_L, a)),
                        ustrip(uconvert(u_L, Λ)),  γ )
