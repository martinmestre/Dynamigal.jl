"""
Customized potential type constructors
"""

"""MilkyWayBovy2014
MWPotential2014 = [
    PowerSphericalPotentialwCutoff(normalize=0.05, alpha=1.8, rc=1.9 / 8.0),
    MiyamotoNagaiPotential(a=3.0 / 8.0, b=0.28 / 8.0, normalize=0.6),
    NFWPotentia
    l(a=2.0, normalize=0.35),
]
https://github.com/jobovy/galpy/blob/b7c3bf055880d21f4b250981acfdd5f0b4f5db09/galpy/potential/mwpotentials.py#L24

Actual values take from Gala code:
https://github.com/adrn/gala/blob/619d167985bdf5091757fa75bb11ba938f5d97fa/src/gala/potential/potential/builtin/special.py#L307
"""
function MilkyWayBovy2014()
    bulge = PowerLawCutoff(m=4501365375.06545*u"Msun", α=1.8, c=1.0*u"kpc")
    disk = MiyamotoNagaiDisk(m=68193902782.346756*u"Msun" , a=3.0*u"kpc", b=280.0*u"pc")
    halo = NFW(m=4.3683325e11*u"Msun", a=16*u"kpc")
    return CompositePotential(bulge,disk,halo)
end
function MilkyWayBovy2014(fac::T) where {T<:Real}
    bulge = PowerLawCutoff(m=4501365375.06545*u"Msun", α=1.8, c=1.0*u"kpc")
    disk = MiyamotoNagaiDisk(m=68193902782.346756*u"Msun" , a=3.0*u"kpc", b=280.0*u"pc")
    halo = NFW(m=fac*4.3683325e11*u"Msun", a=16*u"kpc")
    return CompositePotential(bulge,disk,halo)
end

function MilkyWayPriceWhelan2017()
    nucleus = Hernquist(m=1.71e9*u"Msun", a=0.07*u"kpc")
    bulge = Hernquist(m=5.0e9*u"Msun", a=1.0*u"kpc")
    disk = MiyamotoNagaiDisk(m=6.8e10*u"Msun" , a=3.0*u"kpc", b=280.0*u"pc")
    halo = NFW(m=5.4e11*u"Msun", a=15.62*u"kpc")
    return CompositePotential(nucleus,bulge,disk,halo)
end

""" MilkyWayMosquera2026 customized potential
Defined by Mercedes Mosquera for a paper on SBI orbits in the local group."""
function MilkyWayMosquera2026()
    nucleus = Hernquist(m=1.71e9*u"Msun", a=0.07*u"kpc")
    bulge = Hernquist(m=5.0e9*u"Msun", a=1.0*u"kpc")
    disk = MiyamotoNagaiDisk(m=6.8e10*u"Msun" , a=3.0*u"kpc", b=280.0*u"pc")
    halo = NFW(m=7.41e11*u"Msun", a=17.04*u"kpc")
    return CompositePotential(nucleus,bulge,disk,halo)
end