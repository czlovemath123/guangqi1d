from math import *
import math
import sympy as sp
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import RegularGridInterpolator
from phy_const import *

x_ratio=0.74
y_ratio=0.26

def get_positive_quartic_roots(a, b, c, d, e, threshold=1e-10):
    x = sp.symbols('x')
    eq = a*x**4 + b*x**3 + c*x**2 + d*x + e
    solutions = sp.solve(eq, x)
    
    positive_roots = []
    for sol in solutions:
        num_val = complex(sp.N(sol, 20))
        
        # Your important relative magnitude check
        if (abs(num_val.imag)/abs(num_val.real)) < threshold and num_val.real>0:
            positive_roots.append(float(num_val.real))

    return sorted(positive_roots, reverse=False) if positive_roots else None

def get_positive_cubic_roots(a, b, c, d, threshold=1e-10):
    x = sp.symbols('x')
    eq = a*x**3 + b*x**2 + c*x + d
    solutions = sp.solve(eq, x)

    positive_roots = []
    for sol in solutions:
        # Evaluate numerically with higher precision
        num_val = complex(sp.N(sol, 20))  # 20 digits precision

        # Check if real (imaginary part negligible) and positive
        if (abs(num_val.imag)/abs(num_val.real)) < threshold and num_val.real > 0:
            positive_roots.append(float(num_val.real))

    return sorted(positive_roots, reverse=False) if positive_roots else None

def eradeg(rho,t,x=0.74):
    #calculate erad to eg ratio
    erad=a_rad*t**4
    eg=e_internal(rho,t,x)
    value=erad/eg
    return value

def solve_for_t(rho,eratio,x=0.74,t_low=1e1,t_high=1e6):
    #use bisection method to solve the temperature with given rho and erad to eg ratio
    equation = lambda t: eradeg(rho, t) - eratio
    t_solution = bisect(equation, t_low, t_high)
    return t_solution

def solve_for_t_gamma(rho,eratio,gamma,mu=1.0):
    n=rho/(mu*mh)
    x=eratio*n*kb/(a_rad*(gamma-1))
    t=pow(x,0.3333)
    return t

def zh2(t):
    value=ztr(t,mh2)
    return value
def zh(t):
    value=ztr(t,mh)*exp(-dish/2/kb/t)
    return value
def zhion(t):
    value=ztr(t,mhion)*exp(-(dish+2*ionh)/2/kb/t)
    return value
def ze(t):
    value=ztr(t,me)
    return value
def ztr(t,m):
    value=pow(2*pi*m*kb*t,1.5)/pow(h,3)
    return value
def zheI(t):
    value=ztr(t,mhe)
    return value
def zheII(t):
    value=ztr(t,mhe)*exp(-ionhe1/kb/t)
    return value
def zheIII(t):
    value=ztr(t,mhe)*exp(-(ionhe1+ionhe2)/kb/t)
    return value
def e_internal(rho,t,x):
    epsilon_h2=1.5*kb*t
    epsilon_hi=1.5*kb*t+dish/2
    epsilon_hii=1.5*kb*t+dish/2+ionh
    epsilon_elec=1.5*kb*t
    epsilon_hei=1.5*kb*t
    epsilon_heii=1.5*kb*t+ionhe1
    epsilon_heiii=1.5*kb*t+ionhe1+ionhe2
    nh2,nhi,nhii,nhelec,nhei,nheii,nheiii,nheelec=solve_species_n(rho,t,x_ratio)
    value=nh2*epsilon_h2+nhi*epsilon_hi+nhii*epsilon_hii+nhelec*epsilon_elec+nhei*epsilon_hei   \
        +nheii*epsilon_heii+nheiii*epsilon_heiii+nheelec*epsilon_elec
    return value
def p_internal(rho,t,x):
    nh2,nhi,nhii,nhelec,nhei,nheii,nheiii,nheelec=solve_species_n(rho,t,x_ratio)
    value=(nh2+nhi+nhii+nhelec+nhei+nheii+nheiii+nheelec)*kb*t
    return value

def load_gamma_table(filename):
    # Load the entire file as a 2D numpy array
    raw_data = np.loadtxt(filename)
    
    # First row is log10 density
    rho_array = raw_data[0, :]
    # Second row is log10 temperature
    t_array = raw_data[1, :]
    # The rest is the gamma table
    # Each row in gamma_table corresponds to a temperature in t_array
    # So the shape should be (len(t_array), len(rho_array))
    gamma_table = raw_data[2:, :]
    
    return rho_array, t_array, gamma_table

def interpolate_gamma(logrho, logt, rho_axis, t_axis, table):
    """
    Performs 2D interpolation on the gamma table.
    Note: table structure is (len(t_axis), len(rho_axis))
    """
    # Create the interpolator. 
    # We pass (t_axis, rho_axis) because rows=temp, columns=density
    interp = RegularGridInterpolator((t_axis, rho_axis), table, 
                                     bounds_error=False, fill_value=None)
    
    # The input to the interpolator must be a point (logt, logrho)
    point = np.array([logt, logrho])
    return interp(point)[0]

def adiabatic_cs(rho,t,rho_array, t_array, gamma_table):
    logrho=log10(rho)
    logt=log10(t)
    gamma=interpolate_gamma(logrho,logt,rho_array,t_array,gamma_table)
    p=p_internal(rho,t,x_ratio)
    cs=sqrt(gamma*p/rho)
    return cs

def solve_species_n(rho,t,x):
    #solve the species number of h and he mixture
    y=1.0-x
    rhox=x*rho
    rhoy=y*rho
    th1=10**2.7
    if (rhox>1e-10):
        th2=10**(3.2+0.032*10)
    else:
        th2=10**(3.2+0.032*(log10(rhox)+20))
    if (rhox>1e-10):
        th3=min(10**(3.2+0.038*10+1.25*(log10(rhox)+10)),1.1e6)
    else:
        th3=10**(3.2+0.038*(log10(rhox)+20))
    th4=min(10**(4+0.2*(log10(rhox)+20)),1.1e6)
    the1=10**3.5
    if (rhoy<1e-15):
        the2=10**(4.05+0.02*(log10(rhoy)+20))
    else:
        the2=10**(4+0.03*(log10(rhoy)+20))
    if (rhoy>=1e-5):
        the3=min(10**(5+0.8*log10(rhoy)),1.1e6)
    elif (rhoy>1e-10 and rhoy<1e-5):
        the3=min(10**(4.6+0.08*(log10(rhoy)+10)),1.1e6)
    else:
        the3=10**(4.2+0.04*(log10(rhoy)+20))
    the4=min(10**(4.5+0.2*(log10(rhoy)+20)),1.1e6)
    if (t<=th1):
        nh2,nhi,nhii,nhelec=h_species_state1(rhox)
    elif (t>th4):
        nh2,nhi,nhii,nhelec=h_species_state5(rhox)
    elif (t>th1 and t<=th2):
        nh2,nhi,nhii,nhelec=h_species_state2(rhox,t)
    elif (t>th2 and t<=th3):
        nh2,nhi,nhii,nhelec=h_species_state4(rhox,t)
    else:
        nh2,nhi,nhii,nhelec=h_species_state3(rhox,t)
    if (t<=the1):
        nhei,nheii,nheiii,nheelec=he_species_state1(rhoy)
    elif (t>the4):
        nhei,nheii,nheiii,nheelec=he_species_state5(rhoy)
    elif (t>the1 and t<the2):
        nhei,nheii,nheiii,nheelec=he_species_state2(rhoy,t)
    elif (t>the3 and t<the4):
        nhei,nheii,nheiii,nheelec=he_species_state4(rhoy,t)
    else:
        nhei,nheii,nheiii,nheelec=he_species_state3(rhoy,t)
    return (nh2,nhi,nhii,nhelec,nhei,nheii,nheiii,nheelec)

def h_species_state1(rho):
    #purely molecular
    nh2=rho/mh2
    nhI=0
    nhII=0
    nelec=0
    return (nh2,nhI,nhII,nelec)
def h_species_state2(rho,t):
    #h2 and HI
    qdis=zh(t)*zh(t)/zh2(t)
    nhtot=rho/mh
    nhI=2*qdis*nhtot/(qdis+sqrt(qdis*qdis+8*qdis*nhtot))
    nh2=nhI*nhI/qdis
    nhII=0
    nelec=0
    return (nh2,nhI,nhII,nelec)
def h_species_state3(rho,t):
    #all present
    q1=zhion(t)*ze(t)/zh(t)
    q2=zh(t)*zh(t)/zh2(t)
    a=2.0/q1/q1/q2
    b=0
    c=1.0/q1
    d=1.0
    e=-rho/mh
    roots=get_positive_quartic_roots(a, b, c, d, e, threshold=1e-10)
    value=roots[0]
    nhII=value
    nelec=value
    nhI=nhII**2/q1
    nh2=nhI**2/q2
    return (nh2,nhI,nhII,nelec)
def h_species_state4(rho,t):
    #HI and HII
    qion=zhion(t)*ze(t)/zh(t)
    nhtot=rho/mh
    nhII=2*qion*nhtot/(qion+sqrt(qion*qion+4*qion*nhtot))
    nh2=0
    nhI=nhII*nhII/qion
    nelec=nhII
    return (nh2,nhI,nhII,nelec)
def h_species_state5(rho):
    #fully ionized
    nh2=0
    nhI=0
    nhII=rho/mh
    nelec=rho/mh
    return (nh2,nhI,nhII,nelec)
def he_species_state1(rho):
    #pure HeI
    nheI=rho/mhe
    nheII=0
    nheIII=0
    nheelec=0
    return (nheI,nheII,nheIII,nheelec)
def he_species_state2(rho,t):
    #HeI and HeII
    qheII=zheII(t)*ze(t)/zheI(t)
    nhetot=rho/mhe
    nheII=2*qheII*nhetot/(qheII+sqrt(qheII*qheII+4*qheII*nhetot))
    nheIII=0
    nheI=nheII*nheII/qheII
    nheelec=nheII
    return (nheI,nheII,nheIII,nheelec)
def he_species_state3(rho,t):
    #all present
    q1=zheII(t)*ze(t)/zheI(t)
    q2=zheIII(t)*ze(t)/zheII(t)
    a=q1-4.0*q2
    b=q1*rho/mhe
    c=-q1*q2*(q2+2*rho/mhe)
    d=q1*q2*q2*rho/mhe
    roots=get_positive_cubic_roots(a, b, c, d, threshold=1e-10)
    value=roots[0]
    nheIII=value
    nheelec=2*nheIII/(1-nheIII/q2)
    nheII=nheIII*nheelec/q2
    nheI=nheII*nheelec/q1
    return (nheI,nheII,nheIII,nheelec)
def he_species_state4(rho,t):
    #heII and heIII
    qheIII=zheIII(t)*ze(t)/zheII(t)
    nhetot=rho/mhe
    nheIII=2*nhetot*qheIII/(nhetot+qheIII+sqrt(nhetot**2+qheIII**2+6*qheIII*nhetot))
    nheI=0
    nheelec=nheIII+nhetot
    nheII=nheIII*nheelec/qheIII
    return (nheI,nheII,nheIII,nheelec)
def he_species_state5(rho):
    #pure HeIII
    nheI=0
    nheII=0
    nheIII=rho/mhe
    nheelec=2*nheIII
    return (nheI,nheII,nheIII,nheelec)
