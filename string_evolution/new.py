from __future__ import division
import numpy as np

f=open("COLVAR","r")
lines=f.readlines()
phi=[]
psi=[]
for x in lines:
    phi.append(x.split(' ')[2])
    psi.append(x.split(' ')[3])
f.close()

phi=np.array(phi, dtype=float)
psi=np.array(psi, dtype=float)

phi_new=np.mean(phi)
psi_new=np.mean(psi)

with open("../../phi-psi-ite-2", "a") as myfile:
    myfile.write("%.11f \t"%phi_new)
    myfile.write("%.11f \n"%psi_new)
    
myfile.close()
