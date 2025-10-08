
"""ODE for the Newtonian case: test particle in a potential."""
function _ode(u::AbstractArray{<:Real}, p::P, t::T) where {P<:AbstractPotential, T<:Real}
    return SA[u[siss]..., acceleration(p, u[sis], t)...]
end
# ode(below) is a little faster than _ode(above)
function ode(u::AbstractArray{<:Real}, p::P, t::T) where {P<:AbstractPotential, T<:Real}
    a = acceleration(p, u[sis], t)
    return SVector{6}(u[4], u[5], u[6], a[1], a[2], a[3])
end

"""
ODE for the Newtonian case: system of macro particles.
This function is called when calling <evolve> while using SystemTrait: GenSys.

Benchmark with 20 macroparticles with compound potentials:

-Using acceleration! ( test against ode!)
    Range (min … max):  92.292 s … 98.913 s  ┊ GC (min … max): 4.34% … 4.40%
    Time  (median):     94.966 s             ┊ GC (median):    4.29%
    Time  (mean ± σ):   95.189 s ±  1.897 s  ┊ GC (mean ± σ):  4.31% ± 0.04%
    Memory estimate: 35.77 GiB, allocs estimate: 900456046.

-Using acceleration! ( test against acceleration_c!)
    Range (min … max):  105.302 s … 110.355 s  ┊ GC (min … max): 3.97% … 4.16%
    Time  (median):     107.286 s              ┊ GC (median):    4.08%
    Time  (mean ± σ):   107.894 s ±   1.845 s  ┊ GC (mean ± σ):  4.08% ± 0.07%
    Memory estimate: 35.69 GiB, allocs estimate: 898334683.

-Using acceleration_c! (test against acceleration!)
    Range (min … max):  108.558 s … 113.131 s  ┊ GC (min … max): 4.75% … 4.88%
    Time  (median):     110.515 s              ┊ GC (median):    4.81%
    Time  (mean ± σ):   110.635 s ±   1.223 s  ┊ GC (mean ± σ):  4.81% ± 0.04%
    Memory estimate: 42.26 GiB, allocs estimate: 944679357.

    Conclusion (tested upto 20 macro particles): better to use "ode" with "acceleration!"
        but not so much difference if using "acceleration_c"

The function below is mutation only in the system's acceleration, but no in du argument
"""
function ode(u::AbstractArray{L}, p::MacroParticleSystem, t::T) where {L<:Real,T<:Real}
    n = Integer(length(u)/2)
    return SA[u[n+1:end]..., acceleration!(p, u[begin:n], t)...]
end


"""
ODE for the Newtonian case: system of macro particles.
This function is called when calling <evolve> while using SystemTrait: GenSysMutODE.

Benchmark with 20 macroparticles with compound potentials:
    Test against ode calling acceleration!
    Range (min … max):  107.892 s …  110.167 s  ┊ GC (min … max): 6.10% … 6.08%
    Time  (median):     108.223 s               ┊ GC (median):    6.12%
    Time  (mean ± σ):   108.691 s ± 880.760 ms  ┊ GC (mean ± σ):  6.10% ± 0.02%
    Memory estimate: 56.26 GiB, allocs estimate: 1451240161.

    Conclusion (tested upto 20 macro particles): better to use "ode" above with "acceleration!"
"""
function ode!(du::AbstractArray{L}, u::AbstractArray{L}, p::MacroParticleSystem, t::T) where {L<:Real,T<:Real}
    n = length(p)
    @inbounds for i in 1:n
        # posiciones y velocidades
        xi = @view u[selec(i):selec(i)+2]
        vi = @view u[n*3+selec(i):n*3+selec(i)+2]

        # acumulador de aceleración
        acc_i = zero(MVector{3,L})

        @inbounds for j in 1:n
            if j != i
                xj = @view u[selec(j):selec(j)+2]
                acc_i += acceleration(p[j].pot, xi - xj, t)
            end
        end

        # escribir resultado en du
        du[selec(i):selec(i)+2] .= vi          # derivada de posición = velocidad
        du[n*3+selec(i):n*3+selec(i)+2] .= acc_i  # derivada de velocidad = aceleración
    end
    return nothing
end

"""ODE for the Newtonian case: LargeCloudMW system."""
function ode(u::AbstractArray{L}, p::LargeCloudMW, t::T) where {L<:Real,T<:Real}
    n = Integer(length(u)/2)
    return SA[u[n+1:end]..., acceleration(p, u, t)...]
end