import numpy as np
import matplotlib.pyplot as plt

myfile= open("phe-data.dat", "a")

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
    
    #print(phi)
    for i in range(1,(len(phi)-1)):
        #print(phi[i])
        #myfile.write('/mnt/Data2/avijeet/autoencoder/fresh/string/'+'md'+str(i+1)+'/iteration-'+str(count-1)+'/COLVAR\t'+str(phi[i])+'\t'+str(psi[i])+'\t80\t80\n')
        #myfile.write('/mnt/Data2/avijeet/autoencoder/fresh/string/md'+str(i+1)+'/iteration-'+str(count-1)+'/COLVAR\t'+f"{phi[i]}\t{psi[i]}\t80\t80\n")
        myfile.write(
        f"/media/avijeet/Elements/autoencoder/phenylanine/replica_exchange/md{i+1}/iteration-{count-1}/COLVAR\t"
        f"{phi[i]:.6f}\t{psi[i]:.6f}\t100\t100\t272\n")
    myfile.write('\n')
myfile.close()
