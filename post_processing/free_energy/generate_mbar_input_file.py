import numpy as np
import matplotlib.pyplot as plt

myfile= open("sample_data.dat", "a")

#-- Here, I am reading the simulation performed during the string evolution and using them as input for MBAR
for count in range(31,51):
    phi=[]; psi=[]; 
    f=open("phi-psi-ite-" + str(count), 'r')
    lines=f.readlines()
    for x in lines:
        if not x.strip():
            continue
        phi.append(x.split('\t')[0])
        psi.append(x.split('\t')[1])
    f.close()
    
    phi = np.array(phi, dtype=float)
    psi = np.array(psi, dtype=float)
    
    # Here, add the location of your COLVAR files, followed by the CVs locations, and force constant. Finally, the simulation temperature
    for i in range(1,(len(phi)-1)):
        myfile.write(
        f"/media/avijeet/Elements/autoencoder/phenylanine/replica_exchange/md{i+1}/iteration-{count-1}/COLVAR\t"
        f"{phi[i]:.6f}\t{psi[i]:.6f}\t100\t100\t272\n")
    myfile.write('\n')
myfile.close()
