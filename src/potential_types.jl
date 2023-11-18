"""Potential types"""

@with_kw struct PointMass{T} <: AbstractPotential
    m::T
end
PointMass(m::T) where {T<:Unitful.Mass} = PointMass( ustrip(uconvert(u_M, m)) )


@with_kw struct Plummer{T,D} <: AbstractPotential
    m::T
    b::D
end
Plummer(m::T, b::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
Plummer( ustrip(uconvert(u_M, m)),  ustrip(uconvert(u_L, b)) )


@with_kw struct MiyamotoNagaiDisk{T,D,F} <: AbstractDiskPotential
    m::T
    a::D
    b::F
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
MiyamotoNagaiDisk( ustrip(uconvert(u_M, m)),  ustrip(uconvert(u_L, a)), ustrip(uconvert(u_L, b)) )


@with_kw struct AllenSantillanHalo{M<:U.Mass, L<:U.Length, D<:Real, T<:Real} <: AbstractDiskPotential
    m_u::M
    a_u::L
    Λ_u::L
    γ::D
    m::T
    a::T
    Λ::T
    function AllenSantillanHalo{M,L,D}(m_u, a_u, Λ_u, γ) where {M,L,D}
        m = uconvert(u_M, m_u).val
        a = uconvert(u_L, a_u).val
        Λ = uconvert(u_L, Λ_u).val
        T=typeof(m)
        return new{M,L,D,T}(m_u, a_u, Λ_u, γ, m, a, Λ)
    end
end
AllenSantillanHalo(m_u::M, a_u::L, Λ_u::L, γ::D) where {M,L,D} = AllenSantillanHalo{M,L,D}(m_u, a_u, Λ_u, γ)
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Real,D<:Real,F<:Real,G<:Real} = AllenSantillanHalo(m*u_M, a*u_L, Λ*u_L, γ)