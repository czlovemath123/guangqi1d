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
var_list=[]
for i in range(n_var):
    line=f_var.readline()
    var_list.append(line.rstrip())
print(f"output_var_info.dat: \n{var_list}")
line=f_var.readline()
f_var.close()
filehead=line.rstrip()
nml=f90nml.read('global.data')
xmin=nml['meshinfo']['n_domain'][0]
xmax=nml['meshinfo']['n_domain'][1]
nml2=f90nml.read('problem.data')
mp=nml2['rhd_quantities']['m_planet']
mp=mp*mjupiter
xmin_rj=xmin/rjupiter
xmax_rj=xmax/rjupiter
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
Erad=glb_cell_var(filename,'Erad')
#Erad_int=glb_cell_var(filename,'Erad_int')
Frad=glb_inter_var(filename,'Fradx')
entropy=glb_cell_var(filename,'entropy')
entropy=entropy
#h2=glb_cell_var(filename,'H2')
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
for i in range(nx-1):
    tau_planck[nx-i-2]=tau_planck[nx-i-1]+(x_inter[nx-i]-x_inter[nx-i-1])*(kp[nx-i-2]*rho[nx-i-2]+kp[nx-i-1]*rho[nx-i-1])/2
for i in range(nx-1):
    tau_rosseland[nx-i-2]=tau_rosseland[nx-i-1]+(x_inter[nx-i]-x_inter[nx-i-1])*(kr[nx-i-2]*rho[nx-i-2]+kr[nx-i-1]*rho[nx-i-1])/2
for i in range(nx):
    frad_cell[i]=Frad[i]+(Frad[i+1]-Frad[i])/(x_inter[i+1]-x_inter[i])*(x_center[i]-x_inter[i])
    frad[i]=frad_cell[i]/Erad[i]/c_light
    mfp[i]=1.0/(kr[i]*rho[i])/rjupiter
for i in range(nx+1):
    Frad[i]=Frad[i]*4*pi*x_inter[i]*x_inter[i]/lsun

x_center=x_center/rjupiter
x_inter=x_inter/rjupiter
xmin_jupiter=xmin/rjupiter
xmax_jupiter=xmax/rjupiter
ps_idx=find_post_idx(x_center,temp,erad_temp)
r_ps=x_center[ps_idx]
l_ps=Frad[ps_idx]

dlum=(Frad[-1]-l_ps)*lsun
egv_ps=egv[ps_idx]
rho_ps=rho[ps_idx]
dphi=mp*G/r_ps/rjupiter
mdot=-4*pi*(x_center[-1]*rjupiter)**2*rho[-1]*v[-1]
lacc=mdot*dphi
epsilon_ps=egv_ps/rho_ps
p_ps=p[ps_idx]
pre_idx=ps_idx+10
print(v[pre_idx])
print(0.5*rho[pre_idx]*v[pre_idx]**3*4*pi*(x_center[pre_idx]*rjupiter)**2)
print(lacc,dlum,epsilon_ps*mdot,p_ps/rho_ps)

xticks=[1.5,2,4,20]

fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True,squeeze=True)
print(f'x:\n{x_center}')
print(f'velocity:\n{v/1e5}')
ln1=axes[0,0].plot(x_center,v/1e5,'r-',label=r'$v_{r}$')
axes[0,0].set_xlim(xmin_jupiter,xmax_jupiter)
axes[0,0].set_xscale('log')
axes[0,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,0].set_xlabel(r'$r\ [R_{J}]$')
axes[0,0].set_ylabel(r'$v\ [km\cdot s^{-1}]$')

ln1=axes[1,0].plot(x_center,kr,'r-',label=r'$\kappa_{R}$')
ln2=axes[1,0].plot(x_center,kp,'b-',label=r'$\kappa_{P}$')
lns=ln1+ln2
labs = [l.get_label() for l in lns]
print('luminosity in solar unit: ',l_ps)
print('bottom pressure in bar: ',p[0]/1e6)
axes[1,0].legend(lns,labs,loc=0)
axes[1,0].set_xlim(xmin_jupiter,xmax_jupiter)
axes[1,0].set_xlabel(r'$r\ [R_{J}]$')
axes[1,0].set_xscale('log')
axes[1,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,0].set_ylabel(r'$\kappa\ [cm^{2}\cdot g^{-1}]$')
axes[1,0].set_yscale('log')
axes[1,0].set_ylim(2e-3,2e1)

ln1=axes[2,0].plot(x_inter,1e4*Frad,'k-',label=r'$L_{\rm{r}}$')
axes[2,0].plot(r_ps,l_ps*1e4,'ko',markersize=8,fillstyle='none')
axes[2,0].set_xlim(xmin_jupiter,xmax_jupiter)
axes[2,0].set_xlabel(r'$r\ [R_{J}]$')
axes[2,0].set_xscale('log')
axes[2,0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,0].set_ylabel(r'$L\ [\times10^{-4}\ L_{\odot}]$')

axes[0,1].plot(x_center,rho,'k-',label=r'$\rho$')
axes[0,1].set_xlim(xmin_jupiter,xmax_jupiter)
axes[0,1].set_xscale('log')
axes[0,1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,1].set_xlabel(r'$r\ [R_{J}]$')
axes[0,1].set_ylabel(r'$\rho\ [g\cdot cm^{-3}]$')
axes[0,1].set_yscale('log')


ln1=axes[1,1].plot(x_center,temp,'r-',label=r'$T_{\rm gas}$')
ln2=axes[1,1].plot(x_center,erad_temp,'b-',label=r'$T_{\rm rad}$')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(xmin_jupiter,xmax_jupiter)
axes[1,1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,1].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[1,1].set_ylabel(r'$T\ [K]$')
axes[1,1].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1,1].set_yscale('log')


axes[2,1].plot(x_center,entropy,'k-')
axes[2,1].set_xscale('log')
axes[2,1].set_xlim(xmin_jupiter,xmax_jupiter)
axes[2,1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,1].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[2,1].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2,1].set_ylabel(r'$s\ [k_{B}\cdot m_{H}^{-1}]$')


ln1=axes[0,2].plot(x_center,frad,'r',label=r'$f_{\rm{r}}$')
axes[0,2].set_xscale('log')
axes[0,2].set_xticks(xticks)
axes[0,2].set_xlim(xmin_jupiter,xmax_jupiter)
axes[0,2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,2].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[0,2].set_ylabel(r'$f_{\rm{r}}$')
axes[0,2].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0,2].tick_params(axis='x',which='major',length=5)

fig.tight_layout()
plt.subplots_adjust(hspace=.0)
plt.show()
plt.savefig(cwd+'/pictures/perfectgas'+filenumber+'.png',bbox_inches='tight',dpi=300)
