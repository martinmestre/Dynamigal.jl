@with_kw struct ParticleEnsemble{T<:Real} <: AbstractEnsemble
    m::Vector{T}
    snapshot::Snapshot
end
ParticleEnsemble(m::Vector{M}, t::T, x::Matrix{D}, v::Matrix{F}) where {M<:Real,T<:Real,D<:Real,F<:Real} =
    ParticleEnsemble(m, Snapshot(t, x, v))
ParticleEnsemble(m::Vector{M}, x::Matrix{D}, v::Matrix{F}) where {M<:Real, D<:Real,F<:Real} =
    ParticleEnsemble(m, Snapshot(x, v))
ParticleEnsemble(m::T, t::T, x::Vector{D}, v::Vector{F}) where {T<:Real, D<:Real,F<:Real} =
    ParticleEnsemble(m*ones(size(x)[1]),t, x, v)
ParticleEnsemble(m::T, x::Vector{D}, v::Vector{F}) where {T<:Real, D<:Real,F<:Real} =
    ParticleEnsemble(m*ones(size(x)[1]), x, v)
ParticleEnsemble(x::Vector{D}, v::Vector{F}) where {D<:Real,F<:Real} =
    Particle(1.0ones(size(x)[1]), x, v)
ParticleEnsemble(m::Vector{M}, t::T, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, T<:Unitful.Time,   D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(ustrip(uconvert(𝕦.m, m)), Snapshot(t, x, v))
ParticleEnsemble(m::Vector{M}, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(ustrip(uconvert(𝕦.m, m)), Snapshot(x, v))
ParticleEnsemble(m::M, t::T, x::Vector{D}, v::Vector{F}) where {M<:Unitful.Mass, T<:Unitful.Time, D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(ustrip(uconvert(𝕦.m, m)), Snapshot(t, x, v))
ParticleEnsemble(x::Vector{D}, v::Vector{F}) where {D<:Unitful.Length, F<:Unitful.Velocity} =
    ParticleEnsemble(1.0ones(size(x)[1]), Snapshot(x, v))


