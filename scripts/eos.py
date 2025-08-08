from math import *
from phy_const import *

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
    y=1.0-x
    epsilon_h2=1.5*kb*t
    epsilon_hi=1.5*kb*t+dish/2
    epsilon_hii=1.5*kb*t+dish/2+ionh
    epsilon_elec=1.5*kb*t
    epsilon_hei=1.5*kb*t
    epsilon_heii=1.5*kb*t+ionhe1
    epsilon_heiii=1.5*kb*t+ionhe1+ionhe2
    th1=1000
    th2=pow(10,3.2+0.036*(log10(rho)+20))
    th3=50000
    the1=10**3.5
    if (rho*y<1e-15):
        the2=10**(4.05+0.02*(log10(rho)+20))
    else:
        the2=10**(4+0.03*(log10(rho)+20))
    if (rho>=1e-5):
        the3=min(10**(9+0.8*log10(rho)),1.1e6)
    elif (rho>1e-10 and rho<1e-5):
        the3=min(10**(4.6+0.08*(log10(rho)+10)),1.1e6)
    else:
        the3=10**(4.2+0.04*(log10(rho)+20))
    the4=min(10**(4.5+0.2*(log10(rho)+20)),1.1e6)
    if (t<=th1):
        nh2,nhi,nhii,nhelec=h_species_state1(rho,x)
    elif (t>th1 and t<=th2):
        nh2,nhi,nhii,nhelec=h_species_state2(rho,t,x)
    elif (t>th2 and t<=th3):
        nh2,nhi,nhii,nhelec=h_species_state3(rho,t,x)
    else:
        nh2,nhi,nhii,nhelec=h_species_state4(rho,x)
    if (t<=the1):
        nhei,nheii,nheiii,nheelec=he_species_state1(rho,y)
    elif (t>the4):
        nhei,nheii,nheiii,nheelec=he_species_state4(rho,y)
    elif (t>the1 and t<the2):
        nhei,nheii,nheiii,nheelec=he_species_state2(rho,t,y)
    else:
        nhei,nheii,nheiii,nheelec=he_species_state3(rho,t,y)
    value=nh2*epsilon_h2+nhi*epsilon_hi+nhii*epsilon_hii+nhelec*epsilon_elec+nhei*epsilon_hei   \
        +nheii*epsilon_heii+nheiii*epsilon_heiii+nheelec*epsilon_elec
    value=value
    return value
def p_internal(rho,t,x):
    y=1.0-x
    th1=1000
    th2=pow(10,3.2+0.036*(log10(rho)+20))
    th3=50000
    the1=10**3.5
    if (rho*y<1e-15):
        the2=10**(4.05+0.02*(log10(rho)+20))
    else:
        the2=10**(4+0.03*(log10(rho)+20))
    if (rho>=1e-5):
        the3=min(10**(9+0.8*log10(rho)),1.1e6)
    elif (rho>1e-10 and rho<1e-5):
        the3=min(10**(4.6+0.08*(log10(rho)+10)),1.1e6)
    else:
        the3=10**(4.2+0.04*(log10(rho)+20))
    the4=min(10**(4.5+0.2*(log10(rho)+20)),1.1e6)
    if (t<=th1):
        nh2,nhi,nhii,nhelec=h_species_state1(rho,x)
    elif (t>th1 and t<=th2):
        nh2,nhi,nhii,nhelec=h_species_state2(rho,t,x)
    elif (t>th2 and t<=th3):
        nh2,nhi,nhii,nhelec=h_species_state3(rho,t,x)
    else:
        nh2,nhi,nhii,nhelec=h_species_state4(rho,x)
    if (t<=the1):
        nhei,nheii,nheiii,nheelec=he_species_state1(rho,y)
    elif (t>the4):
        nhei,nheii,nheiii,nheelec=he_species_state4(rho,y)
    elif (t>the1 and t<the2):
        nhei,nheii,nheiii,nheelec=he_species_state2(rho,t,y)
    else:
        nhei,nheii,nheiii,nheelec=he_species_state3(rho,t,y)
    value=(nh2+nhi+nhii+nhelec+nhei+nheii+nheiii+nheelec)*kb*t
    return value
def h_species_state1(rho,x):
    nh2=x*rho/mh2
    nhI=0
    nhII=0
    nelec=0
    return (nh2,nhI,nhII,nelec)
def h_species_state2(rho,t,x):
    qdis=zh(t)*zh(t)/zh2(t)
    nhtot=x*rho/mh
    nhI=2*qdis*nhtot/(qdis+sqrt(qdis*qdis+8*qdis*nhtot))
    nh2=nhI*nhI/qdis
    nhII=0
    nelec=0
    return (nh2,nhI,nhII,nelec)
def h_species_state3(rho,t,x):
    qion=zhion(t)*ze(t)/zh(t)
    nhtot=x*rho/mh
    nhII=2*qion*nhtot/(qion+sqrt(qion*qion+4*qion*nhtot))
    nh2=0
    nhI=nhII*nhII/qion
    nelec=nhII
    return (nh2,nhI,nhII,nelec)
def h_species_state4(rho,x):
    nh2=0
    nhI=0
    nhII=x*rho/mh
    nelec=x*rho/mh
    return (nh2,nhI,nhII,nelec)
def he_species_state1(rho,y):
    nheI=y*rho/mhe
    nheII=0
    nheIII=0
    nheelec=0
    return (nheI,nheII,nheIII,nheelec)
def he_species_state2(rho,t,y):
    qheII=zheII(t)*ze(t)/zheI(t)
    nhetot=y*rho/mhe
    nheII=2*qheII*nhetot/(qheII+sqrt(qheII*qheII+4*qheII*nhetot))
    nheIII=0
    nheI=nheII*nheII/qheII
    nheelec=nheII
    return (nheI,nheII,nheIII,nheelec)
def he_species_state3(rho,t,y):
    q1=zheII(t)*ze(t)/zheI(t)
    q2=zheIII(t)*ze(t)/zheII(t)
    p=np.zeros(4)
    p[0]=(q1-4*q2)
    p[1]=q1*y*rho/mhe
    p[2]=-q1*q2*(q2+y*rho/mhe)
    p[3]=q1*q2*q2*y*rho/mhe
    roots=np.roots(p)
    for i in range(3):
        if (np.isreal(roots[i]) and roots[i]>0):
            value=roots[i]
    nheIII=np.real(value)
    nheelec=2*nheIII/(1-nheIII/q2)
    nheII=nheIII*nheelec/q2
    nheI=nheII*nheelec/q1
    return (nheI,nheII,nheIII,nheelec)
def he_species_state4(rho,y):
    nheI=0
    nheII=0
    nheIII=y*rho/mhe
    nheelec=2*nheIII
    return (nheI,nheII,nheIII,nheelec)
