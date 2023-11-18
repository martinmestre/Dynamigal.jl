"""Constants"""

"""Solver units"""
const lu = UnitsConfig()
const G = ustrip( uconvert(lu.l*lu.v^2/lu.m, u"G") )