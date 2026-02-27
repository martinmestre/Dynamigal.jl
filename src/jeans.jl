"""Velocity dispersions"""

function velocity_dispersion(pot::P; β::R=0.0, r_min::L=0.1, r_max::L=0.0, n_nodes::I=100) where {P<:AbstractPotential, R<:Real, L<:Real, I<:Integer}
    if r_max == 0.0
        g(x) = density(pot,x) - 𝕛.ϵ_ρ
        D(f)= x->gradient(y->f(y),x)[1]
        r_max = find_zero((g, D(g)),  [1.0e-6,1.e3], Roots.Brent())
    end
    @show r_min r_max
    rₐ =range(r_min, r_max, n_nodes)
    σ² = Vector{typeof(rₐ[begin])}(undef, n_nodes)
    vec(r) = [1,1,1]*r/sqrt(3)
    ν(r) = density(pot, vec(r))
    dϕ(r) = norm(acceleration(pot,vec(r)))
    h(u,p) = u^(2β)*ν(u)*dϕ(u)
    for i ∈ eachindex(rₐ)
        prob = IntegralProblem(h, (rₐ[i], Inf))
        σ²[i] = (1/(rₐ[i]^(2β)*ν(rₐ[i])))*solve(prob, QuadGKJL()).u
    end
    itp = interpolate(σ², BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, rₐ)
    return sitp
end