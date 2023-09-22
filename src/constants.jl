"""Constants"""

const u_M = u"Msun"
const u_L = u"kpc"
const u_V = u"km/s"
const u_A = u"km/s^2"
const u_T = uconvert(u"Gyr", 1u"kpc/(km/s)")
const u_Pot = u_V^2
const Gᵤ = uconvert(u_L*u_V^2/u_M, u"G")
const G = Gᵤ.val

@with_kw struct UnitsConfig
    u_M::U.Mass = 1.0u"Msun"
    u_L::U.Length  = 1.0u"kpc"
    u_V::U.Velocity = 1.0u"km/s"
    u_A::U.Acceleration = 1.0u"km/s^2"
    u_T::U.Time = uconvert(u"Gyr", 1u_L/u_V)
    Gᵤ::typeof(u"G") = uconvert(unit(u_L*u_V^2/u_M), u"G")
    G::AbstractFloat = Gᵤ.val
end
