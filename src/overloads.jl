

"""Potential types"""

"""CompositePotential types"""
Base.length(cp::CompositePotential) = length(cp.potentials)
Base.getindex(cp::CompositePotential, i::Int) = cp.potentials[i]
Base.iterate(cp::CompositePotential, state...) = iterate(cp.potentials, state...)
Base.eachindex(cp::CompositePotential) = eachindex(cp.potentials)
Base.getindex(cp::CompositePotential, r::AbstractRange{Int}) = MacroParticleSystem(cp.potentials[r])
Base.getindex(cp::CompositePotential, v::AbstractVector{Int}) = MacroParticleSystem(cp.potentials[v])


function Base.:+(p1::T, p2::U) where {T <: AbstractPotential, U <: AbstractPotential}
    return CompositePotential((p1, p2))
end

function Base.:+(cp::CompositePotential, p::T) where {T <: AbstractPotential}
    return CompositePotential((cp.potentials..., p))
end

function Base.:+(p::T, cp::CompositePotential) where {T <: AbstractPotential}
    return CompositePotential((p, cp.potentials...))
end

function Base.:+(cp1::CompositePotential, cp2::CompositePotential)
    return CompositePotential((cp1.potentials..., cp2.potentials...))
end

"""MacroParticleSystem types"""
Base.length(mps::MacroParticleSystem) = length(mps.macroparticles)
Base.getindex(mps::MacroParticleSystem, i::Int) = mps.macroparticles[i]
Base.iterate(mps::MacroParticleSystem, state...) = iterate(mps.macroparticles, state...)
Base.eachindex(mps::MacroParticleSystem) = eachindex(mps.macroparticles)
Base.getindex(mps::MacroParticleSystem, r::AbstractRange{Int}) = MacroParticleSystem(mps.macroparticles[r])
Base.getindex(mps::MacroParticleSystem, v::AbstractVector{Int}) = MacroParticleSystem(mps.macroparticles[v])


function Base.:+(p1::T, p2::U) where {T <: AbstractMacroParticle, U <: AbstractMacroParticle}
    return MacroParticleSystem((p1, p2))
end

function Base.:+(mps::MacroParticleSystem, p::T) where {T <: AbstractMacroParticle}
    return MacroParticleSystem((mps.macroparticles..., p))
end

function Base.:+(p::T, mps::MacroParticleSystem) where {T <: AbstractMacroParticle}
    return MacroParticleSystem((p, mps.macroparticles...))
end

function Base.:+(mps1::MacroParticleSystem, mps2::MacroParticleSystem)
    return MacroParticleSystem((mps1.macroparticles..., mps2.macroparticles...))
end


"""Base.show for NFW struct"""
function Base.show(io::IO, pot::NFW)
    printstyled(io, "NFW halo:\n"; color=:cyan, bold=true)

    fields = [
        ("m", pot.m),
        ("a", pot.a),
        ("c", pot.c),
        ("m_v", pot.m_v),
        ("r_v", pot.r_v),
        ("ρ₀", pot.ρ₀),
        ("𝔸", pot.𝔸),
        ("cosmos", pot.cosmos)
    ]

    for (name, value) in fields
        print(io, "  ")
        printstyled(io, name; color=:yellow)
        print(io, " = ")
        printstyled(io, string(value); color=:white)
        print(io, "\n")
    end
end