#!/bin/bash
export OMP_NUM_THREADS=6,1
export OMP_STACKSIZE=4G
ulimit -s unlimited
if [ -f geo.gen ] ; then rm -rf geo.gen ; fi
ase convert -i cif -o gen $1 geo.gen
echo '
Geometry = GenFormat {
 <<< geo.gen
}

Hamiltonian = xTB {
  Method = "GFN2-xTB"
  kPointsAndWeights = SuperCellFolding {
    2   0   0
    0   2   0
    0   0   2
    0.5 0.5 0.5
  }
}
Driver = GeometryOptimization {
  LatticeOpt = Yes
}' > dftb_in.hsd
dftb+ > log 
rm -rf dftb_in.hsd
x=$(grep " total energy " log | tail -n1 | awk '{print $3}' | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | bc -lq)
echo "scale=20; 315777.090000000 * $x" | bc -lq
rm -rf input.gen
ase convert -i gen -o cif geo_end.gen geo_end.cif
exit 0
