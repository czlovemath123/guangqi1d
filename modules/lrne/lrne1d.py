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
kb=1.380658E-16
a_rad=7.56e-15
day=86400
G=6.67259e-8
msun=1.99e33
rsun=6.96e10
year=31536000
x=0.74
n=10000


class model_class:
    def __init__(self,ms,rin,time,v,mdot,alpha,n):
        self.ms=ms
        self.rin=rin
        self.time=time
        self.v=v
        self.mdot=mdot
        self.alpha=alpha
        self.t=np.linspace(0,time,n)
        self.dt=time/n*day

def save_model(model,modeldir,eos='gamma'):
    # Define the file path
    file_path = 'models.csv'

    # Define the header and the new row
    header = ['ms','rin','time','v','mdot','alpha']
    new_model = [model.ms,model.rin,model.time,model.v,model.mdot,model.alpha]

    # Check if the file exists
    if not os.path.isfile(file_path):
        # If the file doesn't exist, create it with the header
        df = pd.DataFrame(columns=header)
        df.to_csv(file_path, index=False)

    # Append the new row to the CSV file
    df = pd.DataFrame([new_model], columns=header)
    column_width=6
    df.to_csv(file_path, mode='a', header=False, index=False)
    print(f"model appended: {new_model}")

    os.makedirs(modeldir,exist_ok=True)
    os.makedirs(modeldir+'/out',exist_ok=True)
    os.makedirs(modeldir+'/pictures',exist_ok=True)
    source_file = 'output_var_info.dat'
    destination = modeldir+'/output_var_info.dat'
    shutil.copy(source_file, destination)
    if (eos=='gamma'):
        source_file = 'guangqi'
    elif (eos=='hhe'):
        source_file = 'guangqi_hhe'
    destination = modeldir+'/'+modeldir
    shutil.copy(source_file, destination)
    print(f"Directory '{modeldir}' created successfully!")
    formatted_string = df.to_string(col_space=6, index=False)
    with open(modeldir+'/formatted_data.txt', 'w') as file:
        file.write(formatted_string)
    rho,v,temp,mdot=[np.zeros(n) for _ in range(4)]
    for i in range(n):
        v[i]=model.v*sqrt(2*model.ms*msun*G/(model.rin*rsun))
        temp[i]=model.alpha*mh*v[i]**2/2/kb
        mdot[i]=model.mdot*msun/year
        rho[i]=mdot[i]/4/pi/(model.rin*rsun)**2/v[i]
    data=[model.t*day,v,rho,temp]
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
    radiation_acc=False
    probnml=f90nml.read('problem.data')
    probnml['parameters_1d']['larad']=radiation_acc
    glbnml=f90nml.read('global.data')
    probnml.write(modeldir+'/problem.data',force=True)
    glbnml.write(modeldir+'/global.data',force=True)

model=model_class(ms=6,rin=10,time=1,v=1.01,mdot=1,alpha=0.18,n=10000)
i=3
model_directory='model'+str(i).zfill(2)
#save_model(model,modeldir=model_directory)
save_model(model,modeldir=model_directory,eos='hhe')
save_problem(modeldir=model_directory)
