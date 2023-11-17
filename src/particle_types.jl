"""Particle types"""


@with_kw struct Particle{M<:U.Mass,L<:U.Length,V<:U.Velocity,T<:Real} <: AbstractMacroParticle
    m_u::M
    x_u::Vector{L}
    v_u::Vector{V}
    m::T
    x::Vector{T}
    v::Vector{T}
    function Particle{M,L,V}(m_u, x_u, v_u) where {M,L,V}
        m = uconvert(u_M, m_u).val
        x = ustrip(uconvert.(u_L, x_u))
        v = ustrip(uconvert.(u_V, v_u))
        T=typeof(m)
        return new{M,L,V,T}(m_u, x_u, v_u, m, x, v)
    end
end
Particle(m_u::M, x_u::Vector{L}, v_u::Vector{V}) where {M,L,V} = Particle{M,L,V}(m_u, x_u, v_u)

@with_kw struct TestParticle{L<:U.Length,V<:U.Velocity,T<:Real} <: AbstractTestParticle
x_u::Vector{L}
v_u::Vector{V}
x::Vector{T}
v::Vector{T}
function TestParticle{L,V}(x_u, v_u) where {L,V}
    x = ustrip(uconvert.(u_L, x_u))
    v = ustrip(uconvert.(u_V, v_u))
    T=typeof(x[begin])
    return new{L,V,T}(x_u, v_u, x, v)
end
end
TestParticle(x_u::Vector{L}, v_u::Vector{V}) where {L,V} = TestParticle{L,V}(x_u, v_u)