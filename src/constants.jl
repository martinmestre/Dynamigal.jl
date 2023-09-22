"""Constants"""

"""Solver units"""
const u_M = u"Msun"
const u_L = u"kpc"
const u_V = u"km/s"
const u_A = u"km/s^2"
const u_T = uconvert(u"Gyr", 1u"kpc/(km/s)")
const u_Pot = u_V^2
const Gᵤ = uconvert(u_L*u_V^2/u_M, u"G")
const G = Gᵤ.val
