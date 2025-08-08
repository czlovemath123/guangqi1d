import numpy as np
import os.path
import scipy.integrate as integrate
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib
from math import *
matplotlib.rcParams.update({'font.size': 24})
cwd=os.getcwd()
os.chdir(cwd+'/../tables/dust')
o_sub=1300
p_sub=1500
oliv50=[]
pyro50=[]
pyro60=[]
pyro70=[]
pyro80=[]
with open('oliv50.txt') as f:
    for line in f:
        line=[float(s) for s in line.split()]
        oliv50.append(line)
oliv50=np.transpose(np.asarray(oliv50))
lambda1=oliv50[0]
log_lambda1=np.log10(lambda1)
lambda1_min=np.min(lambda1)
lambda1_max=np.max(lambda1)
o50_kappa=oliv50[1]
o50=interpolate.interp1d(lambda1,o50_kappa,kind='cubic')
with open('pyro50_0.01.txt') as f:
    for line in f:
        line=[float(s) for s in line.split()]
        pyro50.append(line)
pyro50=np.transpose(np.asarray(pyro50))
p50_kappa=pyro50[1]
p50=interpolate.interp1d(lambda1,p50_kappa,kind='cubic')
with open('pyro60_0.01.txt') as f:
    for line in f:
        line=[float(s) for s in line.split()]
        pyro60.append(line)
pyro60=np.transpose(np.asarray(pyro60))
p60_kappa=pyro60[1]
p60=interpolate.interp1d(lambda1,p60_kappa,kind='cubic')
with open('pyro70_0.01.txt') as f:
    for line in f:
        line=[float(s) for s in line.split()]
        pyro70.append(line)
pyro70=np.transpose(np.asarray(pyro70))
p70_kappa=pyro70[1]
p70=interpolate.interp1d(lambda1,p70_kappa,kind='cubic')
with open('pyro80_0.01.txt') as f:
    for line in f:
        line=[float(s) for s in line.split()]
        pyro80.append(line)
pyro80=np.transpose(np.asarray(pyro80))
p80_kappa=pyro80[1]
p80=interpolate.interp1d(lambda1,p80_kappa,kind='cubic')
#x is in micron (\mu m)
def planck_function(x,T):
    co=6.626*3/1.38/10
    x=x/1e4
    return 1.0/pow(x,5)/(exp(co/x/T)-1.0)
def fun_absorb_o50(x,teff):
    kappa=o50(x)
    a=2*pi*planck_function(x,teff)*kappa
    return a
def fun_emit_o50(x,tsub):
    kappa=o50(x)
    a=4*pi*planck_function(x,tsub)*kappa
    return a
def fun_absorb_p50(x,teff):
    kappa=p50(x)
    a=2*pi*planck_function(x,teff)*kappa
    return a
def fun_emit_p50(x,tsub):
    kappa=p50(x)
    a=4*pi*planck_function(x,tsub)*kappa
    return a
def fun_absorb_p60(x,teff):
    kappa=p50(x)
    a=2*pi*planck_function(x,teff)*kappa
    return a
def fun_emit_p60(x,tsub):
    kappa=p60(x)
    a=4*pi*planck_function(x,tsub)*kappa
    return a
def fun_absorb_p70(x,teff):
    kappa=p70(x)
    a=2*pi*planck_function(x,teff)*kappa
    return a
def fun_emit_p70(x,tsub):
    kappa=p70(x)
    a=4*pi*planck_function(x,tsub)*kappa
    return a
def fun_absorb_p80(x,teff):
    kappa=p80(x)
    a=2*pi*planck_function(x,teff)*kappa
    return a
def fun_emit_p80(x,tsub):
    kappa=p80(x)
    a=4*pi*planck_function(x,tsub)*kappa
    return a
teff=2874
o50_flux=integrate.quad(fun_absorb_o50,lambda1_min,20,args=teff,limit=1000)
o50_emit=integrate.quad(fun_emit_o50,lambda1_min,20,args=o_sub,limit=1000)
r_o50=1/sin(acos(1-o50_emit[0]/o50_flux[0]))
print('oliv50 forms around',round(r_o50,2),'stellar radius')
p50_flux=integrate.quad(fun_absorb_p50,lambda1_min,20,args=teff,limit=1000)
p50_emit=integrate.quad(fun_emit_p50,lambda1_min,20,args=p_sub,limit=1000)
r_p50=1/sin(acos(1-p50_emit[0]/p50_flux[0]))
r_pp50=sqrt(1/(1-pow(1-p50_emit[0]/p50_flux[0],2)))
print('pyro50 forms around',round(r_pp50,2),round(r_p50,2),'stellar radius')
p60_flux=integrate.quad(fun_absorb_p60,lambda1_min,20,args=teff,limit=1000)
p60_emit=integrate.quad(fun_emit_p60,lambda1_min,20,args=p_sub,limit=1000)
r_p60=1/sin(acos(1-p60_emit[0]/p60_flux[0]))
r_pp60=sqrt(1/(1-pow(1-p60_emit[0]/p60_flux[0],2)))
print('pyro60 forms around',round(r_pp60,2),round(r_p60,2),'stellar radius')
p70_flux=integrate.quad(fun_absorb_p70,lambda1_min,20,args=teff,limit=1000)
p70_emit=integrate.quad(fun_emit_p70,lambda1_min,20,args=p_sub,limit=1000)
r_p70=1/sin(acos(1-p70_emit[0]/p70_flux[0]))
r_pp70=sqrt(1/(1-pow(1-p70_emit[0]/p70_flux[0],2)))
print('pyro70 forms around',round(r_p70,2),round(r_pp70,2),'stellar radius')
p80_flux=integrate.quad(fun_absorb_p80,lambda1_min,20,args=teff,limit=1000)
p80_emit=integrate.quad(fun_emit_p80,lambda1_min,20,args=p_sub,limit=1000)
r_p80=1/sin(acos(1-p80_emit[0]/p80_flux[0]))
r_pp80=sqrt(1/(1-pow(1-p80_emit[0]/p80_flux[0],2)))
print('pyro80 forms around',round(r_p80,2),round(r_pp80,2),'stellar radius')
b1=np.zeros(lambda1.shape[0])
b2=np.zeros(lambda1.shape[0])
for i in range(lambda1.shape[0]):
    x=lambda1[i]
    b1[i]=planck_function(x,teff)
    b2[i]=planck_function(x,1500)
fig,ax1=plt.subplots(figsize=(16,8))
ax2=ax1.twinx()
ln1=ax1.plot(lambda1,p50_kappa,'k-',linewidth=2,label=r'$\rm{Mg_{0.5}Fe_{0.5}SiO}_{3}$')
ln2=ax1.plot(lambda1,p60_kappa,'k--',linewidth=2,label=r'$\rm{Mg_{0.6}Fe_{0.4}SiO}_{3}$')
ln3=ax1.plot(lambda1,p70_kappa,'k-.',linewidth=2,label=r'$\rm{Mg_{0.7}Fe_{0.3}SiO}_{3}$')
ln4=ax1.plot(lambda1,p80_kappa,'k:',linewidth=2,label=r'$\rm{Mg_{0.8}Fe_{0.2}SiO}_{3}$')
ln5=ax2.plot(lambda1,b1,'b--',linewidth=1,label=r'$B_{\lambda}(2874{\rm K})$')
#ln6=ax2.plot(lambda1,b2,'r--',linewidth=1,label='1500 K')
lns=ln1+ln2+ln3+ln4+ln5
labs=[l.get_label() for l in lns]
ax1.legend(lns,labs,loc=0,ncol=2)
ax1.axvspan(0.1,2.15,alpha=0.5,color='grey')
ax1.axvspan(2.15,2.4,alpha=0.5,color='orange')
ax1.set_xlabel(r'$\lambda\ \rm{[}\mu\rm{m]}$')
ax1.set_ylabel(r'$\kappa_{\lambda}\ \rm{[cm}^{2}\cdot\rm{g}^{-1}\rm{]}}$')
ax2.set_ylabel(r'$I_{\lambda}\ \rm{ [erg}\cdot\rm{cm}^{-2}\cdot\rm{s}^{-1}\rm{]}$')
ax2.set_yscale('log')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(lambda1_min,lambda1_max)
ax1.set_ylim(0.5,5e5)
ax2.set_ylim(1e7,1e19)
fig.tight_layout()
#plt.show()
plt.savefig('dust.pdf',dpi=300)
plt.close()


#r=np.linspace(1.5,10,num=81)
#point_source=np.linspace(0,0,num=81)
#finite_volume=np.linspace(0,0,num=81)
#for i in range(81):
#    x=r[i]
#    theta=asin(1.0/x)
#    finite_volume[i]=2*(1.0-cos(theta))
#    point_source[i]=1/pow(x,2)
#plt.figure(figsize=(16,8))
#plt.plot(r,point_source,'k-',label=r'$r_{0}^{2}/r^{2}$')
#plt.plot(r,finite_volume,'k--',label=r'$2(1-\cos{\theta_{0}})$')
#plt.legend()
#plt.xlabel(r'$r_{0}/r$')
#plt.xlim(1.5,10)
#plt.tight_layout()
##plt.show()
#plt.savefig('pointsource.pdf',dpi=300)
#plt.close()
