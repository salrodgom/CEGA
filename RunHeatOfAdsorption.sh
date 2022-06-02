#!/bin/bash
function check {
 n_max=12
 n=$(ps aux | grep 'simulate' | sed '/grep/d' | wc -l)
 while [ $(echo "${n}>=${n_max}" | bc -lq) == 1 ] ; do
  sleep 10
  n=$(ps aux | grep 'simulate' | sed '/grep/d' | wc -l)
 done
}
gfortran xyzemsemble2def.f90 -o xyzemsemble2def.exe
for folder in dir_* ; do
    molecule=$(echo $folder | sed 's/dir_//g')
    if [ -f ${folder}/${molecule}_Conformer_Ensemble_weights_and_charges.xyz ] && [ ! -d ${folder}/Heat_TAMOF-1 ] ; then
     echo $folder
     mkdir $folder/Heat_TAMOF-1
     cp xyzemsemble2def.exe $folder/.
     cp LOL_scaled.cif force_field_mixing_rules.def $folder/Heat_TAMOF-1/.
     cd $folder
      ln -s ${molecule}_Conformer_Ensemble_weights_and_charges.xyz ${molecule}_all.xyz
      ./xyzemsemble2def.exe -i ${molecule}_all.xyz ; rm -rf xyzemsemble2def.exe ${molecule}_all.xyz
      mv simulation.input pseudo_atoms.def ${molecule}*def Heat_TAMOF-1/.
      cd Heat_TAMOF-1
       check
       nohup simulate > simulate.log &
       sleep 0.5
      cd .. 
     cd ..
    fi
done
