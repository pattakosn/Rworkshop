#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=12
#PBS -q bigmem
#PBS -o out.12.out
#PBS -e out.12.err

cd /home/patkos/evaluation/050
OUT=out.12

#mpiexec Rscript SarahTaxa2DistDDMatrix.r > ${OUT}.Taxa2Dist
mpiexec Rscript SarahTaxa2DistMPI.r > ${OUT}.distMPI
#mpiexec Rscript SarahTaxonDiveMPI.r > ${OUT}.TaxonDive

exit 0 
