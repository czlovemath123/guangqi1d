import sys
sys.path.insert(1,'../../scripts')
from assemble_1d_data import *
from math import *
import matplotlib
import matplotlib.pyplot as plt
import os.path
plt.rcParams.update({'font.size': 14})
np.set_printoptions(threshold=sys.maxsize)
matplotlib.rc('lines',linewidth=2)

def find_post_idx(r,temp,erad_temp):
    n=len(r)
    for i in range(n):
        ps_idx=i
        if (temp[i]-erad_temp[i]>10):
            break
    return i

def center(x):
    n=len(x)-1
    xc=np.zeros(n)
    for i in range(n):
        x1=x[i+1]
        x2=x[i]
        dx=x2-x1
        xm=(x1+x2)/2
        xc[i]=xm+2*xm*dx**2/(12*xm**2+dx**2)
    return xc


fig_format='png'
f_var=open('output_var_info.dat','r')
n_var=int(f_var.readline())
print(n_var)
var_list=[]
for i in range(n_var):
    line=f_var.readline()
    var_list.append(line.rstrip())
print(var_list)
line=f_var.readline()
f_var.close()
filehead=line.rstrip()
nml=f90nml.read('global.data')
xmin=nml['meshinfo']['n_domain'][0]
xmax=nml['meshinfo']['n_domain'][1]
lengthscale=nml['meshinfo']['lengthscale']
xmin=xmin*lengthscale
xmax=xmax*lengthscale
nml2=f90nml.read('problem.data')
mp=nml2['rhd_quantities']['m_planet']
mp=mp*mjupiter
xmin_rj=xmin/lengthscale
xmax_rj=xmax/lengthscale
cwd=os.getcwd()
os.chdir(cwd+'/out')
if (len(sys.argv)>=2):
    i=int(sys.argv[1])
else:
    i=0
filenumber=str(i).zfill(5)
filename=filehead+filenumber+'.h5'
print(filename)


x_center=glb_x_center(filename)
x_inter=glb_x_inter(filename)
#level=glb_level(filename)
rho=glb_cell_var(filename,'rho')
egv=glb_cell_var(filename,'egv')
v=glb_cell_var(filename,'vx')
p=glb_cell_var(filename,'pres')
egv=glb_cell_var(filename,'egv')
Erad=glb_cell_var(filename,'Erad')
Erad_int=glb_cell_var(filename,'Erad_int')
Frad=glb_inter_var(filename,'Fradx')
Ehydro_flux=glb_inter_var(filename,'hydro_energy_xflux')
mass_flux=glb_inter_var(filename,'mass_xflux')
entropy=glb_cell_var(filename,'entropy')
entropy=entropy/NA/kb
h2=glb_cell_var(filename,'H2')
#hii=glb_cell_var(filename,'HII')
kp=glb_cell_var(filename,'Planck')
kr=glb_cell_var(filename,'Rosseland')
temp=glb_cell_var(filename,'temp')
erad_temp=np.power(Erad/arad,0.25)
nx=len(x_center)
tau_planck=1e-5*np.ones(nx)
tau_rosseland=1e-5*np.ones(nx)
frad_cell=np.zeros(nx)
frad=np.zeros(nx)
mfp=np.zeros(nx)
R=np.zeros(nx+1)
gpotential_flux=np.zeros(nx+1)
potential=np.zeros(nx)
for i in range(nx-1):
    tau_planck[nx-i-2]=tau_planck[nx-i-1]+(x_inter[nx-i]-x_inter[nx-i-1])*(kp[nx-i-2]*rho[nx-i-2]+kp[nx-i-1]*rho[nx-i-1])/2
for i in range(1,nx):
    tau_rosseland[i]=tau_rosseland[i-1]+(x_inter[i]-x_inter[i-1])*(kr[i-1]*rho[i-1]+kr[i]*rho[i])/2
for i in range(nx):
    frad_cell[i]=Frad[i]+(Frad[i+1]-Frad[i])/(x_inter[i+1]-x_inter[i])*(x_center[i]-x_inter[i])
    frad[i]=frad_cell[i]/Erad_int[i]/c_light
    mfp[i]=1.0/(kr[i]*rho[i])/lengthscale
for i in range(nx+1):
    Frad[i]=Frad[i]*4*pi*x_inter[i]*x_inter[i]/lsun
    mass_flux[i]=mass_flux[i]*4*pi*x_inter[i]*x_inter[i]
    Ehydro_flux[i]=Ehydro_flux[i]*4*pi*x_inter[i]*x_inter[i]
    gpotential_flux[i]=-G*mp/x_inter[i]*mass_flux[i]
#print(Frad*lsun)
#print(Ehydro_flux)
#print(gpotential_flux)
#print(mass_flux)
#exit()

x_center=x_center/lengthscale
x_inter=x_inter/lengthscale
x_center_log=np.log10(x_center)
x_inter_log=np.log10(x_inter)
xmin_rp=log10(xmin/lengthscale)
xmin_rp=-0.01
xmax_rp=log10(xmax/lengthscale)
ps_idx=find_post_idx(x_center,temp,erad_temp)
r_ps=x_center_log[ps_idx]
l_ps=Frad[ps_idx]

print(tau_rosseland[ps_idx])

xticks=[1.1,1.7,4,20]
#fig,axes=plt.subplots(1,1,figsize=(12,10),sharex=True,squeeze=True)
#axes.plot(x_center,egv/rho,'k-')
#plt.show()

totalflux=Ehydro_flux+gpotential_flux
fig,axes=plt.subplots(1,1,figsize=(12,10),sharex=True,squeeze=True)
axes.plot(x_inter_log,(Frad-Frad[ps_idx])*lsun,'k-',label='Frad')
#axes.plot(x_inter,Ehydro_flux,'r-')
#axes.plot(x_inter,gpotential_flux,'b-')
axes.plot(x_inter_log,-(totalflux-totalflux[ps_idx]),'k--',label='total')
axes.set_xlim(xmin_rp,xmax_rp)
plt.legend()
#axes.set_xscale('log')
plt.savefig('../flux.png',dpi=100)
plt.show()

#exit()

#fig,axes=plt.subplots(2,2,figsize=(12,10),sharex=True,squeeze=True,gridspec_kw={'height_ratios': [1,1,1,1.5,1]})
fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True,squeeze=True)

#ax2=axes[0,0].twinx()
ln1=axes[0,0].plot(x_center_log,v/1e5,'r-',label=r'$v_{r}$')
#ln2=ax2.plot(x_center,rho,'k-',linewidth=1,label=r'$\rho$')
#lns=ln1+ln2
#labs=[l.get_label() for l in lns]
axes[0,0].set_xlim(xmin_rp,xmax_rp)
#axes[0,0].set_xscale('log')
axes[0,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,0].set_xlabel(r'$r\ [R_{J}]$')
axes[0,0].set_ylabel(r'$v\ [km\cdot s^{-1}]$')
#ax2.set_ylabel(r'$\rho\ [g\cdot cm^{-3}]$')
#ax2.set_yscale('log')
#axes[0,0].legend(lns,labs,loc=0)

#ax2=axes[1,0].twinx()
ln1=axes[1,0].plot(x_center_log,kr,'r-',label=r'$\kappa_{R}$')
ln2=axes[1,0].plot(x_center_log,kp,'b-',label=r'$\kappa_{P}$')
lns=ln1+ln2
labs = [l.get_label() for l in lns]
print('luminosity in solar unit:',l_ps)
print('bottom pressure in bar: ',p[0]/1e6)
axes[1,0].legend(lns,labs,loc=0)
axes[1,0].set_xlim(xmin_rp,xmax_rp)
axes[1,0].set_xlabel(r'$r\ [R_{J}]$')
#axes[1,0].set_xscale('log')
axes[1,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,0].set_ylabel(r'$\kappa\ [cm^{2}\cdot g^{-1}]$')
axes[1,0].set_yscale('log')
axes[1,0].set_ylim(2e-3,2e1)

ln1=axes[2,0].plot(x_inter_log,1e4*Frad,'k-',label=r'$L_{\rm{r}}$')
axes[2,0].plot(r_ps,l_ps*1e4,'ko',markersize=8,fillstyle='none')
axes[2,0].set_xlim(xmin_rp,xmax_rp)
axes[2,0].set_xlabel(r'$r\ [R_{J}]$')
#axes[2,0].set_xscale('log')
axes[2,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,0].set_ylabel(r'$L\ [\times10^{-4}\ L_{\odot}]$')





#ax2=axes[0,1].twinx()
#ln1=axes[2,0].plot(x_center,h2,'r-',linewidth=1,label=r'$\chi_{\rm H_2}$')
#ln2=axes[2,0].plot(x_center,hii,'b-',linewidth=1,label=r'$\chi_{\rm H^{+}}$')
#ln3=ax2.plot(x_center,level,'k-',linewidth=1,label='level')
#lns=ln1+ln2+ln3
#labs=[l.get_label() for l in lns]
#axes[2].axvspan(2.4,4.8,alpha=0.5,color='pink')
#axes[2].text(3.2,0.2,'H2 dissociation')
#axes[2,0].set_xlim(xmin_rp,xmax_rp)
#axes[2,0].set_xscale('log')
#axes[2,0].legend(lns,labs,loc=5)
#axes[2,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#axes[2,0].set_xlabel(r'$r\ [R_{J}]$')
#axes[2,0].set_ylabel(r'$\chi$')
#ax2.set_ylabel('level')


axes[0,1].plot(x_center_log,rho,'k-',label=r'$\rho$')
axes[0,1].set_xlim(xmin_rp,xmax_rp)
#axes[0,1].set_xscale('log')
axes[0,1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,1].set_xlabel(r'$r\ [R_{J}]$')
axes[0,1].set_ylabel(r'$\rho\ [g\cdot cm^{-3}]$')
axes[0,1].set_yscale('log')


#ax2=axes[3].twinx()
ln1=axes[1,1].plot(x_center_log,temp,'r-',label=r'$T_{\rm gas}$')
ln2=axes[1,1].plot(x_center_log,erad_temp,'b-',label=r'$T_{\rm rad}$')
#ln3=ax2.plot(x_center,entropy,'k-',linewidth=1,label=r'$s$')
#lns=ln1+ln2+ln3
#labs=[l.get_label() for l in lns]
#axes[3].axvspan(2.4,4.8,alpha=0.5,color='orange')
#axes[3].text(2.6,1000,'endothermic')
#axes[1,1].set_xscale('log')
axes[1,1].set_xlim(xmin_rp,xmax_rp)
axes[1,1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,1].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[1,1].set_ylabel(r'$T\ [K]$')
axes[1,1].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,1].set_yscale('log')
#axes[0,1].legend(lns,labs,loc=0)
#ax2.set_ylabel(r'$s\ [k_{B}\cdot m_{H}^{-1}]$')


axes[2,1].plot(x_center_log,entropy,'k-')
#axes[2,1].set_xscale('log')
axes[2,1].set_xlim(xmin_rp,xmax_rp)
axes[2,1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,1].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[2,1].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,1].set_ylabel(r'$s\ [k_{B}\cdot m_{H}^{-1}]$')


ln1=axes[0,2].plot(x_center_log,p/1e6,'k',label='pressure')
#axes[0,2].set_xscale('log')
axes[0,2].set_yscale('log')
axes[0,2].set_xticks(xticks)
axes[0,2].set_xlim(xmin_rp,xmax_rp)
axes[0,2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,2].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[0,2].set_ylabel('p [bar]')
axes[0,2].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,2].tick_params(axis='x',which='major',length=5)

ln1=axes[1,2].plot(x_center_log,frad,'r',label=r'$f_{\rm{r}}$')
#axes[1,2].set_xscale('log')
axes[1,2].set_xticks(xticks)
axes[1,2].set_xlim(xmin_rp,xmax_rp)
axes[1,2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,2].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[1,2].set_ylabel(r'$f_{\rm{r}}$')
axes[1,2].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,2].tick_params(axis='x',which='major',length=5)

ln1=axes[2,2].plot(x_center_log,h2,'r',label='h2')
#axes[2,2].set_xscale('log')
axes[2,2].set_xticks(xticks)
axes[2,2].set_xlim(xmin_rp,xmax_rp)
axes[2,2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,2].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[2,2].set_ylabel(r'$\chi_{H2}$')
axes[2,2].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,2].tick_params(axis='x',which='major',length=5)


fig.tight_layout()
##fig.suptitle(r'$M_{p}=M_{J},\dot{M}_{p}=10^{-2}M_{\oplus}\cdot\rm{yr}^{-1}$')
#fig.subplots_adjust(top=0.95)
plt.subplots_adjust(hspace=.0)
#
#
#
##a=plt.axes([0.55,0.37,0.1,0.08])
##a.plot(x_axis,erad_temp,'b-',linewidth=3)
##a.plot(x_axis,temp,'r-',linewidth=1)
##b=a.twinx()
##b.plot(x_axis,entropy,'k-',linewidth=1)
##a.set_xticks([1.85,1.87])
##a.set_yticks([2000,7000])
###b.set_yticks([12,14,16,18,20,22])
##a.set_xlim(1.85,1.87)
##a.set_ylim(2000,7000)
##b.set_ylim(15,25)
#
#
#
plt.savefig(cwd+'/sim'+filenumber+'.png',bbox_inches='tight',dpi=300)
print(cwd+'sim'+filenumber+'.png')
##plt.close()
plt.show()


