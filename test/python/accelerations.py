"""Python functions to test Julia code"""
import numpy as np


g = 4.498502151469552e-6 # kpc Gyr^{-1} Msun^{-1}
# g = 4.300923924e-6  # kpc (km/s)^2 M_sun^{-1}...This unit is not used in GalacticDynamics.jl


class Plummer:
    """Plummer Sphere"""

    def __init__(self, m, b):
        """Init."""
        self.m = m
        self.b = b

    def accel(self, x,y,z):
        fac = -g*self.m/(self.b*self.b +x*x + y*y + z*z )**(1.5)
        return np.array([fac*x, fac*y, fac*z])


class MiyamotoNagai:
    """Miyamoto-Nagai disk."""

    def __init__(self, M, a, b):
        """Init."""
        self.M = M
        self.a = a
        self.b = b

    def accel(self, x, y, z):
        """Acceleration."""
        fac = -g*self.M / (x*x + y*y + (self.a + np.sqrt(self.b**2 + z*z))**2)**(1.5)
        return np.array([fac*x, fac*y, fac*z*(1.0+self.a/np.sqrt(self.b**2+z*z))])


class AllenSantillan:
    """Allen-Santillan halo"""

    def __init__(self, M_h, r_h, L, gam):
        """Init."""
        self.M_h = M_h
        self.r_h = r_h
        self.L = L
        self.gam = gam

    def f(self,x):
        return 1.0+(x/self.r_h)**(self.gam-1)

    def accel(self, x,y,z):
        r=np.sqrt(x*x+y*y+z*z)
        if(isinstance(r, np.float64)):
            if(r<self.L):
                fac = g*self.M_h/self.r_h*(1-1/self.f(r))/r**2
            else:
                fac = g*self.M_h*(self.L/self.r_h)**self.gam/self.f(self.L)/r**3
        else:
            fac = np.zeros(len(r))
            for i in range(0,len(r)):
                if(r[i]<self.L):
                    fac[i] = g*self.M_h/self.r_h*(1-1/self.f(r[i]))/r[i]**2
                else:
                    fac[i] = g*self.M_h*(self.L/self.r_h)**self.gam/self.f(self.L)/r[i]**3
        return np.array([fac*x, fac*y, fac*z])

