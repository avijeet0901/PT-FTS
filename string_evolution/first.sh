

#-----Required input
#-----phi-psi-ite-1, with row for each string image and column of value of CVs
#-----All the intermediate configurations in conf folder
#-----gromacs mdp and topology files

for i in {1..27}
do 
mkdir md$i
cd md$i
cp ../phi-psi-ite-1 ./
gmx_mpi grompp -f ../grompp.mdp -c ../../conf/md$i.gro -p ../topol.top -o md$i.tpr

awk "NR==$i{print;exit}" phi-psi-ite-1 > test.dat
phi=$(awk 'NR == 1 {printf "%f\n",$1}' test.dat)
psi=$(awk 'NR == 1 {printf "%f\n",$2}' test.dat)
rm test.dat

cat >plumed.dat << EOF #Update this section according to types and number of CVs
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

gmx_mpi mdrun -s md$i.tpr -deffnm md$i -plumed plumed.dat -nsteps 1000000 -v

cd ..
done
