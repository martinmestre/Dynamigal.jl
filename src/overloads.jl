

"""Potential types"""
# function Base.show(io::IO, pot::T) where {T <: AbstractPotential}
#     fields = fieldnames(T)
#     if isempty(fields)
#         print(io, T.name.name, "()")
#     else
#         field_str = join(["$f=$(getfield(pot, f))" for f in fields], ", ")
#         print(io, T.name.name, "(", field_str, ")")
#     end
# end

# function Base.show(io::IO, ::MIME"text/plain", pot::T) where {T <: AbstractPotential}
#     # Misma implementación que show(io::IO, pot)
#     fields = fieldnames(T)
#     if isempty(fields)
#         print(io, T.name.name, "()")
#     else
#         field_str = join(["$f=$(getfield(pot, f))" for f in fields], ", ")
#         print(io, T.name.name, "(", field_str, ")")
#     end
# end

"""CompositePotential types"""
Base.length(cp::CompositePotential) = length(cp.potentials)
Base.getindex(cp::CompositePotential, i::Int) = cp.potentials[i]
Base.iterate(cp::CompositePotential, state...) = iterate(cp.potentials, state...)
Base.eachindex(cp::CompositePotential) = eachindex(cp.potentials)

# function Base.show(io::IO, cp::CompositePotential)
#     print(io, "CompositePotential($(length(cp)) components)")
# end

# function Base.show(io::IO, ::MIME"text/plain", cp::CompositePotential)
#     n = length(cp)
#     println(io, "CompositePotential with $n components:")
#     for (i, p) in enumerate(cp)
#         println(io, "  [$i] ", p)
#     end
# end

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

print("estoy aca")