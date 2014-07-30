#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=2:ppn=12
#PBS -q bigmem
#PBS -o out.24.out
#PBS -e out.24.err

cd /home/patkos/evaluation/010
OUT=out.24

mpiexec Rscript SarahTaxa2DistDDMatrix.r > ${OUT}.Taxa2Dist
mpiexec Rscript SarahTaxa2DistMPI.r > ${OUT}.distMPI
mpiexec Rscript SarahTaxonDiveMPI.r > ${OUT}.TaxonDive

exit 0 
