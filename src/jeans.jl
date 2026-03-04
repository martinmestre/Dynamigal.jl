"""Velocity dispersions"""

function velocity_dispersion(pot::P; β::R=0.0, r_min::L=0.1, r_max::T=0.0, n_nodes::I=100) where {P, R<:Real, L<:Real, T<:Real, I<:Integer}
    vec(r) = [1,0,1]*r/sqrt(2)
    if r_max == 0.0
        g(x) = density(pot,vec(x)) - 𝕛.ϵ_ρ
        D(f)= x->gradient(y->f(y),x)[1]
        r_max = find_zero((g, D(g)),  [1.0e-6,1.e3], Roots.Brent())
    end
    @show r_min r_max
    rₐ =range(r_min, r_max, n_nodes)
    σ = Vector{typeof(rₐ[begin])}(undef, n_nodes)
    vec(r) = [1,0,1]*r/sqrt(2)
    ν(r) = density(pot, vec(r))
    dϕ(r) = norm(acceleration(pot,vec(r)))
    h(u,p) = u^(2β)*ν(u)*dϕ(u)
    for i ∈ eachindex(rₐ)
        prob = IntegralProblem(h, (rₐ[i], Inf))
        𝕀 = solve(prob,QuadGKJL()).u
        σ[i] = sqrt( (1/(rₐ[i]^(2β)*ν(rₐ[i])))*𝕀 )
    end
    itp = interpolate(σ, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, rₐ)
    return sitp
end