@with_kw struct ParticleEnsemble{T<:Real} <: AbstractEnsemble
    m::Vector{T}
    snapshot::Snapshot
end
ParticleEnsemble(m::Vector{M}, t::T, x::Matrix{D}, v::Matrix{F}) where {M<:Real,T<:Real,D<:Real,F<:Real} =
    ParticleEnsemble(m, Snapshot(t, x, v))
ParticleEnsemble(m::Vector{M}, x::Matrix{D}, v::Matrix{F}) where {M<:Real, D<:Real,F<:Real} =
    ParticleEnsemble(m, Snapshot(x, v))
ParticleEnsemble(m::Vector{M}, t::T, x::Matrix{D}, v::Matrix{F}) where {M<:Unitful.Mass, T<:Unitful.Time,   D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(ustrip(uconvert(𝕦.m, m)), Snapshot(t, x, v))
ParticleEnsemble(m::Vector{M}, x::Matrix{D}, v::Matrix{F}) where {M<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(ustrip(uconvert(𝕦.m, m)), Snapshot(x, v))
ParticleEnsemble(m::M, t::T, x::Matrix{D}, v::Matrix{F}) where {M<:Unitful.Mass, T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(ustrip(uconvert(𝕦.m, m)), Snapshot(t, x, v))
ParticleEnsemble(x::Matrix{D}, v::Matrix{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(1.0ones(size(x)[1]), Snapshot(x, v))


@with_kw struct TestParticleEnsemble{T<:Real} <: AbstractEnsemble
    snapshot::Snapshot
end
TestParticleEnsemble(t::T, x::Matrix{D}, v::Matrix{F}) where {T<:Real,D<:Real,F<:Real} =
    TestParticleEnsemble(Snapshot(t, x, v))
TestParticleEnsemble(x::Matrix{D}, v::Matrix{F}) where {D<:Real,F<:Real} =
    TestParticleEnsemble(Snapshot(x, v))
TestParticleEnsemble(t::T, x::Matrix{D}, v::Matrix{F}) where {T<:Unitful.Time,   D<:Unitful.Length, F<:Unitful.Velocity} =
    TestParticleEnsemble(Snapshot(t, x, v))
TestParticleEnsemble(x::Matrix{D}, v::Matrix{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    TestParticleEnsemble(Snapshot(x, v))
TestParticleEnsemble(x::Matrix{D}, v::Matrix{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    TestParticleEnsemble(Snapshot(x, v))

