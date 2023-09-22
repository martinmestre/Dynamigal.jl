"""Configuration structs"""

@with_kw struct SolverConfig
    solver::supertype(Vern9) = Vern9()
    abstol::Float64 = 0.5e-8
    reltol::Float64 = 5.0e-8
end