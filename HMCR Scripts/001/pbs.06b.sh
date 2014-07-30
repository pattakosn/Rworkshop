#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=6
#PBS -q bigmem
#PBS -o out.06b.out
#PBS -e out.06b.err

cd /home/patkos/evaluation/001
OUT=out.06b

mpiexec Rscript SarahTaxa2DistDDMatrix.r > ${OUT}.Taxa2Dist
mpiexec Rscript SarahTaxa2DistMPI.r > ${OUT}.distMPI
mpiexec Rscript SarahTaxonDiveMPI.r > ${OUT}.TaxonDive

exit 0 
