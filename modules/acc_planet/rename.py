import sys
import os

p = sys.argv[1:]

mdot = p[0]
mp = p[1]
rp = p[2]
lint = p[3] 
old = "out"
new = f"/home/azha/Documents/guangqi_out_0306/simplified/{mdot}_{mp}_{rp}_{lint}"
os.rename(old,new)




