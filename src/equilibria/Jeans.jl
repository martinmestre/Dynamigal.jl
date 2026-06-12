"""Velocity dispersion from spherical Jeans equation"""
function velocity_dispersion(pot::P, r::T, β::R) where {P,T<:Real, R<:Real}
    vec(s) = [1,1,√2]*s/2
    ν(s) = density(pot, vec(s))
    dϕ(s) = norm(acceleration(pot,vec(s)))
    h(s,p) = s^(2β)*ν(s)*dϕ(s)
    prob = IntegralProblem(h, (r, Inf))
    𝕀 = solve(prob, QuadGKJL()).u
    return sqrt( (1/(r^(2β)*ν(r)))*𝕀 )
end

function velocity_dispersion(pot::P, r::T, β::F) where {P,T<:Real, F<:Function}
    h_β(s,p) = β(s)/s
    prob_β(s) = IntegralProblem(h_β, (1, s))
    γ(s) = 2*solve(prob_β(s),QuadGKJL()).u

    vec(s) = [1,1,√2]*s/2
    ν(s) = density(pot, vec(s))
    dϕ(s) = norm(acceleration(pot,vec(s)))
    h(s,p) = exp(γ(s))*ν(s)*dϕ(s)
    prob = IntegralProblem(h, (r, Inf))
    𝕀 = solve(prob, QuadGKJL()).u
    return sqrt( (1/(exp(γ(r))*ν(r)))*𝕀 )
end

function velocity_dispersion(pot::P, β::F; r_min::L=0.1, r_max::T=0.0, n_nodes::I=200) where {P, F, L<:Real, T<:Real, I<:Integer}
    vec(u) = [1,1,√2]*u/2
    if r_max == 0.0
        g(x) = density(pot,vec(x)) - 𝕛.ϵ_ρ
        D(f)= x->gradient(y->f(y),x)[1]
        r_max = find_zero((g, D(g)),  [1.0e-6,1.e3], Roots.Brent())
    end
    rₐ =range(r_min, r_max, n_nodes)
    σ = Vector{typeof(rₐ[begin])}(undef, n_nodes)
    for i ∈ eachindex(rₐ)
        σ[i] = velocity_dispersion(pot, rₐ[i], β)
    end
    return cubic_interp(rₐ, σ; bc=ZeroCurvBC(), extrap=NoExtrap())
end

