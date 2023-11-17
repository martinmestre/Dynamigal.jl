"""Potential types"""

@with_kw struct PointMass{M<:U.Mass, T<:Real} <: AbstractPotential
    m_u::M
    m::T
    function PointMass{M}(m_u) where {M}
        m = uconvert(u_M, m_u).val
        T=typeof(m)
        return new{M,T}(m_u, m)
    end
end
PointMass(m_u::M) where {M} = PointMass{M}(m_u)
PointMass(m::T) where {T<:Real} = PointMass(m*u_M)

@with_kw struct Plummer{M<:U.Mass, L<:U.Length, T<:Real} <: AbstractPotential
    m_u::M
    b_u::L
    m::T
    b::T
    function Plummer{M,L}(m_u, b_u) where {M,L}
        m = uconvert(u_M, m_u).val
        b = uconvert(u_L, b_u).val
        T=typeof(m)
        return new{M,L,T}(m_u, b_u, m, b)
    end
end
Plummer(m_u::M, b_u::L) where {M,L} = Plummer{M,L}(m_u, b_u)
Plummer(m::T, b::D) where {T<:Real,D<:Real} = Plummer(m*u_M, b*u_L)

@with_kw struct MiyamotoNagaiDisk{M<:U.Mass, L<:U.Length, T<:Real} <: AbstractDiskPotential
    m_u::M
    a_u::L
    b_u::L
    m::T
    a::T
    b::T
    function MiyamotoNagaiDisk{M,L}(m_u, a_u, b_u) where {M,L}
        m = uconvert(u_M, m_u).val
        a = uconvert(u_L, a_u).val
        b = uconvert(u_L, b_u).val
        T=typeof(m)
        return new{M,L,T}(m_u, a_u, b_u, m, a, b)
    end
end
MiyamotoNagaiDisk(m_u::M, a_u::L, b_u::L) where {M,L} = MiyamotoNagaiDisk{M,L}(m_u, a_u, b_u)
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Real, D<:Real, F<:Real} = Plummer(m*u_M, a*u_L, b*u_L)


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