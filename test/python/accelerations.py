"""Python functions to test Julia code"""

"""Allen-Santillan Acceleration"""
def grad_halo(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    if(isinstance(r, np.float64)):
        if(r<L):
            fac = g*M_h/r_h*(1-1/f(r))/r**2
        else:
            fac = g*M_h*(L/r_h)**gam/f(L)/r**3
    else:
        fac = np.zeros(len(r))
        for i in range(0,len(r)):
            if(r[i]<L):
                fac[i] = g*M_h/r_h*(1-1/f(r[i]))/r[i]**2
            else:
                fac[i] = g*M_h*(L/r_h)**gam/f(L)/r[i]**3
    return np.array([fac*x, fac*y, fac*z])