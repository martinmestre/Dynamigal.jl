"""Potential types"""


@with_kw struct TimeDependent{T<:Real} <: AbstractPotential
    m::T
    @assert m>0 "m must be possitive"
end
TimeDependent(m::T) where {T<:Unitful.Mass} = TimeDependent( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct Kepler{T<:Real} <: AbstractPotential
    m::T
    @assert m>0 "m must be possitive"
end
Kepler(m::T) where {T<:Unitful.Mass} = Kepler( ustrip(uconvert(𝕦.m, m)) )


@with_kw struct Plummer{T<:Real,D<:Real} <: AbstractPotential
    m::T
    a::D
    @assert (m>0 && a>0) "all fields should be possitive"
end
Plummer(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Plummer( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct Hernquist{T<:Real,D<:Real} <: AbstractPotential
    m::T
    a::D
    @assert m>0 && a>0 "all fields should be possitive"
end
Hernquist(m::T, a::D) where {T<:Unitful.Mass, D<:Unitful.Length} =
    Hernquist( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)) )


@with_kw struct MiyamotoNagaiDisk{T<:Real,D<:Real,F<:Real} <: AbstractDiskPotential
    m::T
    a::D
    b::F
    @assert m>0 && a>0 && b>0 "all fields should be possitive"
end
MiyamotoNagaiDisk(m::T, a::D, b::F) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length} =
    MiyamotoNagaiDisk( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)), ustrip(uconvert(𝕦.l, b)) )


@with_kw struct AllenSantillanHalo{T<:Real,D<:Real,F<:Real,G<:Real} <: AbstractHaloPotential
    m::T
    a::D
    Λ::F
    γ::G
    @assert m>0 && a>0 && Λ>0  "fields m, a, Λ should be possitive"
end
AllenSantillanHalo(m::T, a::D, Λ::F, γ::G) where {T<:Unitful.Mass, D<:Unitful.Length, F<:Unitful.Length, G<:Real} =
    AllenSantillanHalo( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)),
                        ustrip(uconvert(𝕦.l, Λ)),  γ )



# NFW

f_nfw(x::T) where {T<:Real} = log(1+x)-x/(1+x)

function r_vir_nfw(m; 𝕔=𝕔)
    ρ = ustrip( uconvert(𝕦.m/𝕦.l^3, 𝕔.ρ_c))
    r = (m/(200*ρ*4.0/3.0*π))^(1.0/3.0)
    return r
end
r_vir_nfw(m::M; 𝕔=𝕔) where {M<:Unitful.Mass} = r_vir_nfw(adimensional(m); 𝕔=𝕔)


@with_kw struct NFW{T<:Real, F<:Real, D<:Real, C<:AbstractConfig} <: AbstractHaloPotential
    @assert m>0 && a>0  "all fields should be possitive"
    m::T  # virial mass: M(r)
    a::F  # scale radius: a=r/c
    𝕔::C = 𝕔
    r::D = r_vir_nfw(m; 𝕔=𝕔) # virial radius
    c::D = r/a # concentration: c=r/a
    𝔸::D = f_nfw(c)
end
NFW(m::T, a::F; 𝕔=𝕔) where {T,F} = NFW(; m=m, a=a, 𝕔=𝕔)
NFW(m::M, a::L) where {M<:Unitful.Mass, L<:Unitful.Length} =
    NFW( ustrip(uconvert(𝕦.m, m)),  ustrip(uconvert(𝕦.l, a)))


function NFW(m::M, c::T; 𝕔=𝕔) where {M<:Unitful.Mass, T<:Real}
    m = adimensional(m)
    @assert m>0 && c>0  "all fields should be possitive"
    r = r_vir_nfw(m; 𝕔=𝕔)  # virial radius
    a = r/c
    𝔸 = f_nfw(c)
    return NFW(m, a, 𝕔, r, c, 𝔸)
end

function concentration(p::NFW; 𝕔=𝕔)
    ρ = ustrip( uconvert(𝕦.m/𝕦.l^3, 𝕔.ρ_c))
    r = (p.m/(200*ρ*4.0/3.0*π))^(1.0/3.0)  # virial radius
    return r/p.a
end
