#!/bin/bash

temps=(272 302 335 372 413 450)

for ite in {2..50} #loop over number of iterations “N”
do

gro=$(( ite-1 ))

python equispacing.py

mv phi-psi-ite-1 ./phi-psi-ite/phi-psi-ite-$ite

mv phi-psi-ite-2 phi-psi-ite-1

for AT in {2..26}
do
mkdir ./md$AT/iteration-$ite/
cp phi-psi-ite-1 ./md$AT/iteration-$ite/.
cd ./md$AT
cp ./iteration-$gro/md$gro.gro ./iteration-$ite/.
cp new.py ./iteration-$ite/.
cd ./iteration-$ite

for i in "${temps[@]}"; do
    gmx_mpi grompp -f ../../${i}.mdp -c md$gro.gro -p ../../topol.top -o ${i}.tpr
done

awk "NR==$AT{print;exit}" phi-psi-ite-1 > test.dat
phi=$(awk 'NR == 1 {printf "%f\n",$1}' test.dat)
psi=$(awk 'NR == 1 {printf "%f\n",$2}' test.dat)
rm test.dat

cat >plumed.dat << EOF
ene: ENERGY
phi: TORSION ATOMS=5,7,9,25 
psi: TORSION ATOMS=7,9,25,27

chi1:  TORSION ATOMS=7,9,11,14
chi2:  TORSION ATOMS=9,11,14,15

restraint-phi: RESTRAINT ARG=phi KAPPA=100 AT=$phi
restraint-psi: RESTRAINT ARG=psi KAPPA=100 AT=$psi

# monitor the two variables and the bias potential from the two restraints
PRINT STRIDE=1 ARG=phi,psi,ene,chi2,chi1,restraint-phi.bias,restraint-psi.bias FILE=COLVAR
EOF

# Create replica directories and copy tpr files
for i in "${!temps[@]}"; do
    mkdir -p replica_${i}
    mv ${temps[$i]}.tpr replica_${i}/topol.tpr
    cp plumed.dat replica_${i}/
done

mpirun -np 6 gmx_mpi mdrun -multidir replica_0 replica_1 replica_2 replica_3 replica_4 replica_5 -deffnm topol -replex 500 -plumed plumed.dat -ntomp 3 -dlb no -v -gpu_id 0 -nb gpu -bonded cpu 

cp ./replica_0/COLVAR.0 ./COLVAR 
cp ./replica_0/topol.gro ./md$ite.gro
sed -i '/^#/d' COLVAR

python new.py

cd ../..
done
done
