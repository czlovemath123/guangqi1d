import numpy as np
from math import *
import os.path
cwd=os.getcwd()
mh=1.6733e-24
kb=1.380658e-16
maw=1.3
m=400
n=601
gamm=1.05
rhomin=1e-15
rhomax=1e-10
tmin=1e2
tmax=2.5e4
f1=open('gamma.txt','w')
f2=open('eg.txt','w')
f3=open('maw.txt','w')
f4=open('p.txt','w')
dlogrho=(log10(rhomax)-log10(rhomin))/(m-1)
dlogt=(log10(tmax)-log10(tmin))/(n-1)
rho=np.asarray(range(m))
rho=log10(rhomin)+dlogrho*rho
t=np.asarray(range(n))
t=log10(tmin)+dlogt*t
for i in range(m):
    f1.write(str('%18.8e'%rho[i]))
f1.write('\n')
for i in range(n):
    f1.write(str('%18.8e'%t[i]))
f1.write('\n')
for i in range(n):
    for j in range(m):
        f1.write(str('%18.8e'%gamm))
    f1.write('\n')
for i in range(m):
    f2.write(str('%18.8e'%rho[i]))
f2.write('\n')
for i in range(n):
    f2.write(str('%18.8e'%t[i]))
f2.write('\n')
for i in range(n):
    for j in range(m):
        eg=(10**rho[j])*(10**t[i])/mh/maw*kb/(gamm-1)
        f2.write(str('%18.8e'%eg))
    f2.write('\n')
for i in range(m):
    f3.write(str('%18.8e'%rho[i]))
f3.write('\n')
for i in range(n):
    f3.write(str('%18.8e'%t[i]))
f3.write('\n')
for i in range(n):
    for j in range(m):
        f3.write(str('%18.8e'%maw))
    f3.write('\n')
for i in range(m):
    f4.write(str('%18.8e'%rho[i]))
f4.write('\n')
for i in range(n):
    f4.write(str('%18.8e'%t[i]))
f4.write('\n')
for i in range(n):
    for j in range(m):
        p=(10**rho[j])*(10**t[i])/mh/maw*kb
        f4.write(str('%18.8e'%p))
    f4.write('\n')
f1.close()
f2.close()
f3.close()
f4.close()
#meta table for perfect gas
nrho=501
npp=701
negm=701
rhomin=1e-15
rhomax=1e-10
tmin=1e2
tmax=2.5e4
pmin=rhomin*kb*tmin/mh/maw
pmax=rhomax*kb*tmax/mh/maw
egmmin=kb*tmin/mh/maw/(gamm-1)
egmmax=kb*tmax/mh/maw/(gamm-1)
dlogrho=(log10(rhomax)-log10(rhomin))/(nrho-1)
dlogp=(log10(pmax)-log10(pmin))/(npp-1)
dlogegm=(log10(egmmax-log10(egmmin)))/(negm-1)
meta_rho=np.asarray(range(nrho))
meta_rho=log10(rhomin)+dlogrho*meta_rho
meta_p=np.asarray(range(npp))
meta_p=log10(pmin)+dlogp*meta_p
meta_egm=np.asarray(range(negm))
meta_egm=log10(egmmin)+dlogegm*meta_egm
f1=open('arhop.txt','w')
f2=open('egmrhop.txt','w')
f3=open('gammarhop.txt','w')
f4=open('prhoegm.txt','w')
for i in range(nrho):
    f1.write(str('%18.8e'%meta_rho[i]))
f1.write('\n')
for i in range(npp):
    f1.write(str('%18.8e'%meta_p[i]))
f1.write('\n')
for i in range(nrho):
    for j in range(npp):
        a=sqrt(gamm*10**meta_p[j]/10**meta_rho[i])
        f1.write(str('%18.8e'%a))
    f1.write('\n')
for i in range(nrho):
    f2.write(str('%18.8e'%meta_rho[i]))
f2.write('\n')
for i in range(npp):
    f2.write(str('%18.8e'%meta_p[i]))
f2.write('\n')
for i in range(nrho):
    for j in range(npp):
        egm=log10(10**meta_p[j]/10**meta_rho[i]/(gamm-1))
        f2.write(str('%18.8e'%egm))
    f2.write('\n')
for i in range(nrho):
    f3.write(str('%18.8e'%meta_rho[i]))
f3.write('\n')
for i in range(npp):
    f3.write(str('%18.8e'%meta_p[i]))
f3.write('\n')
for i in range(nrho):
    for j in range(npp):
        f3.write(str('%18.8e'%gamm))
    f3.write('\n')
for i in range(nrho):
    f4.write(str('%18.8e'%meta_rho[i]))
f4.write('\n')
for i in range(negm):
    f4.write(str('%18.8e'%meta_egm[i]))
f4.write('\n')
for i in range(nrho):
    for j in range(negm):
        p=log10(10**meta_egm[j]*10**meta_rho[i]*(gamm-1))
        f4.write(str('%18.8e'%p))
    f4.write('\n')
f1.close()
f2.close()
f3.close()
f4.close()
