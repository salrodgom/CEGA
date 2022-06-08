#!/bin/bash
# Create XYZ file form SMILES code (add hydrogen atoms): 
# http://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html
export OMP_NUM_THREADS=8,1
export OMP_STACKSIZE=4G
ulimit -s unlimited

export name=$(echo $1 | sed 's/\.xyz//g')
folder=dir_$name
if [ ! -d $folder ] ; then 
mkdir $folder 
cd $folder
mv ../$1 .
echo "0.0" > .CHRG
xtb $1 --charge 0 --opt extreme > Step_00.log                                         # GFN2-xTB
mv xtbopt.xyz $1
crest $1 -chrg 0 -alpb water --entropy > Step_0.log                                   # Remove/Add --entropy to use or not the iMTD-sMTD workflow, faster.
crest -forall crest_conformers.xyz -chrg 0 -alpb water > Step_1.log         
xtb crest_best.xyz --charge 0 --ohess extreme --alpb water > Step_2.log
crest crest_best.xyz -cregen crest_conformers.xyz -chrg 0 -alpb water > Step_3.log    # Sort conformers
echo "
" > top ; cat charges >> top ; paste crest_conformers.xyz.sorted top > ${name}_Conformer_Ensemble_weights_and_charges.xyz
  rm -rf PROP 
  cd ..
 fi
exit 0
