"""Velocity dispersions"""

# function velocity_dispersion(pot::P; β::R=0.0, r_max::R=0.0) where {P<:AbstractPotential, R<:Real}
#     if r_max == 0.0
#         g(x) = density(pot,x) - 𝕛.ρ
#         D(f)= x->gradient(y->f(y),x)[1]
#         r_max = find_zero((g,D(g)),  [1.0e-6,200.0], Roots.Brent())
#     end
#     @show r_max


# end