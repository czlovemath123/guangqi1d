import sys
import os

p = sys.argv[1:]

mdot = p[0]
mp = p[1]
rp = p[2]
lint = p[3] 
old = "out"
root = os.environ.get('GUANGQI_OUT_ROOT',
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), 'simplified'))
os.makedirs(root, exist_ok=True)
new = os.path.join(root, f"{mdot}_{mp}_{rp}_{lint}")
os.rename(old,new)




