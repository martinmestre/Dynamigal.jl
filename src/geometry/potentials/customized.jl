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
"""
# function MilkyWayBovy2014()
#     bulge = PowerLawCutoff(m=, a=)
#     disk = MiyamotoNagaiDisk(m= , a=3.0, b=0.28)
#     halo = NFW(m= , a=)
#     return CompositePotential(bulge,disk,halo)
# end