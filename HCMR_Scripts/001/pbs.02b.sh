#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=2
#PBS -q bigmem
#PBS -o out.02b.out
#PBS -e out.02b.err

cd /home/patkos/evaluation/001
OUT=out.02b

mpiexec Rscript SarahTaxa2DistDDMatrix.r > ${OUT}.Taxa2Dist
mpiexec Rscript SarahTaxa2DistMPI.r > ${OUT}.distMPI
mpiexec Rscript SarahTaxonDiveMPI.r > ${OUT}.TaxonDive

exit 0 
