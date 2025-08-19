from math import *
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import f90nml
import copy
import shutil
plt.rcParams.update({'font.size': 20})

amu=1.66053886e-24
mh2=3.3466E-24
mh=1.6733e-24
me=9.1093897e-28
mhion=mh-me
mhe=6.646481526e-24
ionh=2.18e-11
dish=7.17e-12
ionhe1=3.94e-11
ionhe2=8.72e-11
kb=1.380658E-16
h=6.6260755E-27
a_rad=7.56e-15
day=86400
yr=3.154e7
G=6.67259e-8
msun=1.99e33
rsun=6.96e10
year=31536000
x=0.74
n=10000
r=10*rsun
ms=6*msun

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

def vej(t,mod):
    if (t<mod.t1):
        if (t<mod.x2*mod.t1):
            value=mod.v1_floor1+(mod.v1-mod.v1_floor1)*(t/(mod.x2*mod.t1))
        else:
            value=mod.v1_floor2+(mod.v1-mod.v1_floor2)*exp(-(mod.x1*(t-mod.x2*mod.t1)/mod.t1)**2) 
    else:
        value=mod.v2+(mod.v1_floor2+(mod.v1-mod.v1_floor2)*exp(-(mod.x1*(1.0-mod.x2))**2)-mod.v2)*exp(-(t-mod.t1)/mod.dt2)-mod.v2decline*(t-mod.t1)/mod.t2
    value=value*sqrt(2*mod.ms*G/r)
    return value

def massloss(t,mod):
    dt=mod.t1/20
    dm=mod.mdot2-mod.mdot1-mod.mdot1rise
    if (t>mod.t1):
        x=(mod.t1-t)/dt
        value=mod.mdot2-dm*exp(x)/(exp(x)+1)-mod.mdot2decline*(t-mod.t1)/mod.t2
    else:
        x=(t-mod.t1)/dt
        value=mod.mdot1+dm*exp(x)/(exp(x)+1)+mod.mdot1rise*(t/mod.t1)**2
    return value

class model_class:
    def __init__(self,ms,t1,t2,dt2,v1,v1_floor1,v1_floor2,v2,x1,x2,mdot1,mdot2,mdot1rise,mdot2decline,v2decline,alpha,n):
        self.ms=ms
        self.t1=t1
        self.t2=t2
        self.dt2=dt2
        self.v1=v1
        self.v1_floor1=v1_floor1
        self.v1_floor2=v1_floor2
        self.v2=v2
        self.x1=x1
        self.x2=x2
        self.mdot1=mdot1
        self.mdot2=mdot2
        self.mdot1rise=mdot1rise
        self.mdot2decline=mdot2decline
        self.v2decline=v2decline
        self.alpha=alpha
        self.t=np.linspace(0,t1+t2,n)
        self.dt=(t1+t2)/n*day

def plot_figure(models):
    m=len(models)
    v=np.zeros((n,m))
    vkin=np.zeros((n,m))
    temp=np.zeros((n,m))
    rhox=np.zeros((n,m))
    mdot=np.zeros((n,m))
    culmdot=np.zeros((n,m))
    eg=np.zeros((n,m))
    erad=np.zeros((n,m))
    eradeg=np.zeros((n,m))
    for j in range(m):
        for i in range(n):
            v[i,j]=vej(models[j].t[i],models[j])
            temp[i,j]=models[j].alpha*mh*v[i,j]**2/2/kb
            mdot[i,j]=massloss(models[j].t[i],models[j])
            if (i==0):
                culmdot[i,j]=mdot[i,j]*models[j].dt/yr
            else:
                culmdot[i,j]=culmdot[i-1,j]+mdot[i,j]*models[j].dt/yr
            rho=mdot[i,j]*msun/year/4/pi/r**2/v[i,j]
            rhox[i,j]=rho
            eg[i,j]=e_internal(rho,temp[i,j],x)
            erad[i,j]=a_rad*temp[i,j]**4
            vkin[i,j]=sqrt(2*(eg[i,j]+erad[i,j])/rho)
            eradeg[i,j]=erad[i,j]/eg[i,j]

    rej0=models[0].ms*G/max(vkin[:,0])**2/rsun
    print('max source kinetic energy in kms:',max(vkin[:,0]/1e5))
    print('max ejecta launch radii in rsun:',rej0)

    vesc=sqrt(2*ms*G/r)
    fig,axes=plt.subplots(3,2,figsize=(16,8),sharex=True,squeeze=True)
    axes[0,0].plot(models[0].t,v[:,0]/1e5,'k-')
    axes[0,0].hlines(vesc/1e5,xmin=0,xmax=30,color='red',linestyle='solid',linewidth=2,label=r'$v_{\rm{esc}}$')
    axes[0,0].legend(fontsize=20)
    axes[0,0].set_ylim(410,510)
    axes[0,0].set_ylabel(r'$v_{\rm{ej}}$'+' [km/s]')
    axes[0,0].axvspan(0,1.3,alpha=0.5,color='orange')
    axes[1,0].plot(models[0].t,vkin[:,0]/1e5,'k-')
    axes[1,0].axvspan(0,1.3,alpha=0.5,color='orange')
    tmin=0
    tmax=15
    axes[1,0].set_ylim(110,450)
    axes[1,0].set_ylabel(r'$v_{\rm{source}}$'+' [km/s]')
    axes[2,0].plot(models[0].t,temp[:,0]/1e5,'k-',label=r'$T_{\rm{outflow}}$')
    axes[2,0].axvspan(0,1.3,alpha=0.5,color='orange')
    axes[2,0].set_ylim(1.6,2.7)
    axes[2,0].set_xlabel('t [day]')
    axes[2,0].set_xlim(tmin,tmax)
    axes[2,0].set_ylabel(r'$T\ [\times10^{5}K]$')
    axes[0,1].plot(models[0].t,mdot[:,0],'k-',label=r'$\dot{m}_{\rm{outflow}}$')
    axes[0,1].set_ylabel(r'$\dot{M}\ [M_{\odot}\cdot\rm{yr}^{-1}]$')
    axes[0,1].axvspan(0,1.3,alpha=0.5,color='orange')
    axes[1,1].plot(models[0].t,culmdot[:,0],'k-',label=r'$\dot{m}_{\rm{outflow}}$')
    axes[1,1].axvspan(0,1.3,alpha=0.5,color='orange')
    axes[1,1].set_ylim(-0.001,0.051)
    axes[1,1].set_ylabel(r'$\Delta M\ [M_{\odot}]$')
    axes[2,1].plot(models[0].t,eradeg[:,0],'k-',label=r'$\mathcal{E}/e_{\rm{g}}$')
    axes[2,1].axvspan(0,1.3,alpha=0.5,color='orange')
    axes[2,1].set_ylabel(r'$\mathcal{E}/e_{\rm{g}}$')
    axes[2,1].set_yscale('log')
    axes[2,1].set_xlabel('t [day]')
    axes[2,1].set_xlim(tmin,tmax)
    axes[2,1].set_ylim(0.2,15)
    plt.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.savefig('boundcond.png',bbox_inches='tight',dpi=300)
    #plt.savefig('boundcond.pdf',bbox_inches='tight',dpi=300)
    plt.show()

def create_model(t1=1.3,t2=8,dt2=1.6,v1=1.06,v1_floor1=1,v1_floor2=0.94,v2=0.87,x1=6,
    x2=0.55,mdot1=0.16,mdot2=1.3,mdot1rise=0.2,mdot2decline=0,v2decline=-0.01,alpha=0.18,n=10000):
    #t1=peak ejecta duration, t2=plateau ejecta duration
    model=model_class(ms,t1,t2,dt2,v1,v1_floor1,v1_floor2,v2,x1,x2,mdot1,mdot2,mdot1rise,mdot2decline,v2decline,alpha,n)
    return model

def save_model(model,modeldir):
    # Define the file path
    file_path = 'models.csv'

    # Define the header and the new row
    header = ['t1','t2','dt2','v1','v1_floor1','v1_floor2','v2','x1','x2','mdot1',
        'mdot2','mdot1rise','mdot2decline','v2decline','alpha']
    new_model = [model.t1,model.t2,model.dt2,model.v1,model.v1_floor1,model.v1_floor2,model.v2,
        model.x1,model.x2,model.mdot1,model.mdot2,model.mdot1rise,model.mdot2decline,model.v2decline,model.alpha]

    # Check if the file exists
    df = pd.DataFrame(columns=header)

    # Append the new row to the CSV file
    df = pd.DataFrame([new_model], columns=header)
    column_width=18
    df.to_csv(file_path, mode='a', header=False, index=False)
    print(f"model appended: {new_model}")

    #if (modeldir==''):
    #    print('need to give a modeldir')
    #    exit()
    os.makedirs(modeldir,exist_ok=True)
    os.makedirs(modeldir+'/out',exist_ok=True)
    source_file = 'output_var_info.dat'
    destination = modeldir+'/output_var_info.dat'
    shutil.copy(source_file, destination)
    source_file = 'guangqi'
    destination = modeldir+'/'+modeldir
    shutil.copy(source_file, destination)
    print(f"Directory '{modeldir}' created successfully!")
    formatted_string = df.to_string(col_space=12, index=False)
    with open(modeldir+'/formatted_data.txt', 'w') as file:
        file.write(formatted_string)
    v=np.zeros(n)
    vkin=np.zeros(n)
    temp=np.zeros(n)
    rhox=np.zeros(n)
    mdot=np.zeros(n)
    culmdot=np.zeros(n)
    eg=np.zeros(n)
    erad=np.zeros(n)
    eradeg=np.zeros(n)

    for i in range(n):
        v[i]=vej(model.t[i],model)
        temp[i]=model.alpha*mh*v[i]**2/2/kb
        mdot[i]=massloss(model.t[i],model)
        if (i==0):
            culmdot[i]=mdot[i]*model.dt/yr
        else:
            culmdot[i]=culmdot[i-1]+mdot[i]*model.dt/yr
        rho=mdot[i]*msun/year/4/pi/r**2/v[i]
        rhox[i]=rho
        eg[i]=e_internal(rho,temp[i],x)
        erad[i]=a_rad*temp[i]**4
        vkin[i]=sqrt(2*(eg[i]+erad[i])/rho)
        eradeg[i]=erad[i]/eg[i]

    data=[model.t*day,v,rhox,temp]
    data=np.transpose(data)
    mm=np.shape(data)[0]
    nn=np.shape(data)[1]
    headers=['time','v','rho','temp']
    with open(modeldir+'/bcinput.dat', 'w') as file:
        for header in headers:
            file.write(f"{header:>16}")  # A16 format (16-character width, left-aligned)
        file.write("\n")  # Newline after headers
        for i in range(mm):
            for j in range(nn):
                # Write each number in ES16.8E2 format
                file.write(f"{data[i, j]:16.8E}")
            file.write("\n")  # Newline after each row

def save_problem(modeldir):
    probnml=f90nml.read('problem.data')
    glbnml=f90nml.read('global.data')
    probnml['parameters_1d']['m_star']=ms/msun
    probnml['parameters_1d']['lfld_mom']=True
    probnml['parameters_1d']['larad']=True
    probnml['parameters_1d']['record_length']=2
    probnml['parameters_1d']['petsc_qratio']=1.2
    probnml['parameters_1d']['petsc_iter']=6
    probnml['parameters_1d']['floor_tauR']=1e-2
    probnml['parameters_1d']['rho_floor']=1e-17
    probnml.write(modeldir+'/problem.data',force=True)
    glbnml['meshinfo']['n_domain']=[10,4000,0,0]
    glbnml['meshinfo']['lengthscale']=6.96e10
    glbnml['meshinfo']['tfinal']=50
    glbnml['meshinfo']['timescale']=86400
    glbnml['meshinfo']['CFL']=0.95
    glbnml['meshinfo']['nrefine_region']=1
    glbnml['phyinfo']['llnx']=True
    glbnml['phyinfo']['xgeo_h']=25
    glbnml['phyinfo']['llny']=False
    glbnml['phyinfo']['lam_con']=False
    glbnml['global_parameters']['nd']=1
    glbnml['global_parameters']['nx']=512
    glbnml['global_parameters']['ny']=1
    glbnml['global_parameters']['blk_size_nx']=64
    glbnml['global_parameters']['blk_size_ny']=1
    glbnml['refinement'][0]['refine_xmin']=0
    glbnml['refinement'][0]['refine_xmax']=10.1
    glbnml['refinement'][0]['refine_ymin']=0
    glbnml['refinement'][0]['refine_ymax']=0
    glbnml['refinement'][0]['level']=4
    glbnml.write(modeldir+'/global.data',force=True)

def model_suite(paras,i):
    model=create_model(t1=paras[0],t2=paras[1],dt2=paras[2],v1=paras[3],v1_floor1=paras[4],
        v1_floor2=paras[5],v2=paras[6],x1=paras[7],x2=paras[8],mdot1=paras[9],mdot2=paras[10],
        mdot1rise=paras[11],mdot2decline=paras[12],v2decline=paras[13],alpha=paras[14],n=paras[15])
    i=int(i)
    model_directory='model'+str(i).zfill(2)
    save_model(model,modeldir=model_directory)
    save_problem(modeldir=model_directory)
    return model

model=model_class(ms=ms,t1=1.1,t2=10,dt2=1.6,v1=0.99,v1_floor1=0.97,v1_floor2=0.96,v2=0.92,x1=6,x2=0.65,mdot1=0.16,mdot2=1,mdot1rise=0.2,mdot2decline=0.3,v2decline=0.04,alpha=0.018,n=10000)
i=99
model_directory='model'+str(i).zfill(2)
save_model(model,modeldir=model_directory)
save_problem(modeldir=model_directory)
modelss=[model]
plot_figure(modelss)
