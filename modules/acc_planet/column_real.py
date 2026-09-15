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
x_center=glb_cell_var(filename,'x_center')
x_inter=glb_inter_var(filename,'mesh_x')
level=glb_level(filename)
rho=glb_cell_var(filename,'rho')
v=glb_cell_var(filename,'vx')
p=glb_cell_var(filename,'pres')
Erad=glb_cell_var(filename,'Erad')
Frad=glb_inter_var(filename,'Fradx')
entropy=glb_cell_var(filename,'entropy')
entropy=entropy/NA/kb
h2=glb_cell_var(filename,'H2')
hii=glb_cell_var(filename,'HII')
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
    frad_cell[i]=(Frad[i]+Frad[i+1])/2
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

fig,axes=plt.subplots(5,1,figsize=(8,16),sharex=True,squeeze=True,gridspec_kw={'height_ratios': [1,1,1,1.5,1]})

ax2=axes[0].twinx()
ln1=axes[0].plot(x_center,v/1e5,'r-',linewidth=1,label=r'$v_{r}$')
ln2=ax2.plot(x_center,rho,'k-',linewidth=1,label=r'$\rho$')
lns=ln1+ln2
labs=[l.get_label() for l in lns]
axes[0].set_xlim(xmin_jupiter,xmax_jupiter)
axes[0].set_xscale('log')
axes[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[0].set_xlabel(r'$r\ [R_{J}]$')
axes[0].set_ylabel(r'$v\ [km\cdot s^{-1}]$')
ax2.set_ylabel(r'$\rho\ [g\cdot cm^{-3}]$')
ax2.set_yscale('log')
axes[0].legend(lns,labs,loc=0)

ax2=axes[1].twinx()
ln1=axes[1].plot(x_center,kr,'r-',linewidth=1,label=r'$\kappa_{R}$')
ln2=axes[1].plot(x_center,kp,'b-',linewidth=1,label=r'$\kappa_{P}$')
ln3=ax2.plot(x_inter,1e4*Frad,'k-',linewidth=1,label=r'$L_{\rm{r}}$')
lns=ln1+ln2+ln3
labs = [l.get_label() for l in lns]
print('luminosity in solar unit:',l_ps)
ax2.plot(r_ps,l_ps*1e4,'ko',markersize=8,fillstyle='none')
axes[1].legend(lns,labs,loc=0)
axes[1].set_xlim(xmin_jupiter,xmax_jupiter)
axes[1].set_xlabel(r'$r\ [R_{J}]$')
axes[1].set_xscale('log')
axes[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[1].set_ylabel(r'$\kappa\ [cm^{2}\cdot g^{-1}]$')
axes[1].set_yscale('log')
axes[1].set_ylim(2e-3,2e1)
ax2.set_ylabel(r'$L\ [\times10^{-4}\ L_{\odot}]$')

ax2=axes[2].twinx()
ln1=axes[2].plot(x_center,h2,'r-',linewidth=1,label=r'$\chi_{\rm H_2}$')
ln2=axes[2].plot(x_center,hii,'b-',linewidth=1,label=r'$\chi_{\rm H^{+}}$')
ln3=ax2.plot(x_center,level,'k-',linewidth=1,label='level')
lns=ln1+ln2+ln3
labs=[l.get_label() for l in lns]
#axes[2].axvspan(2.4,4.8,alpha=0.5,color='pink')
#axes[2].text(3.2,0.2,'H2 dissociation')
axes[2].set_xlim(xmin_jupiter,xmax_jupiter)
axes[2].set_xscale('log')
axes[2].legend(lns,labs,loc=5)
axes[2].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[2].set_xlabel(r'$r\ [R_{J}]$')
axes[2].set_ylabel(r'$\chi$')
ax2.set_ylabel('level')

ax2=axes[3].twinx()
ln1=axes[3].plot(x_center,temp,'r-',linewidth=1,label=r'$T_{\rm gas}$')
ln2=axes[3].plot(x_center,erad_temp,'b-',linewidth=1,label=r'$T_{\rm rad}$')
ln3=ax2.plot(x_center,entropy,'k-',linewidth=1,label=r'$s$')
lns=ln1+ln2+ln3
labs=[l.get_label() for l in lns]
#axes[3].axvspan(2.4,4.8,alpha=0.5,color='orange')
#axes[3].text(2.6,1000,'endothermic')
axes[3].set_xscale('log')
axes[3].set_xlim(xmin_jupiter,xmax_jupiter)
axes[3].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[3].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[3].set_ylabel(r'$T\ [K]$')
axes[3].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[3].set_yscale('log')
axes[3].legend(lns,labs,loc=0)
ax2.set_ylabel(r'$s\ [k_{B}\cdot m_{H}^{-1}]$')


ax2=axes[4].twinx()
ln1=axes[4].plot(x_center,frad,'r',linewidth=1,label=r'$f_{\rm{r}}$')
labs=[l.get_label() for l in lns]
#axes[4].axvspan(1.7,1.86,alpha=0.5,color='grey')
#axes[4].text(1.7,0.2,'radiative zone',rotation=90)
axes[4].set_xscale('log')
axes[4].set_xticks([1.1,1.7,2,3,4,20])
axes[4].legend(lns,labs,loc=0)
axes[4].set_xlim(xmin_jupiter,xmax_jupiter)
axes[4].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[4].set_xlabel(r'$r\ [R_{\rm{J}}]$')
axes[4].set_ylabel(r'$f_{\rm{r}}$')
axes[4].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes[4].tick_params(axis='x',which='major',length=5)
ax2.set_ylabel(r'$\lambda_{\rm{mfp}}/r$')
ax2.set_yscale('log')
ax2.set_yticks([1e-3,1,1e3])
ax2.set_ylim([1e-4,4e3])


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
plt.savefig(cwd+'/realgas'+filenumber+'.png',bbox_inches='tight',dpi=300)
##plt.close()
plt.show()
