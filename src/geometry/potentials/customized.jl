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