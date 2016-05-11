#!/bin/bash

#SBATCH -J oofdtd_strong
#SBATCH -o oofdtd_strong.o%j
#SBATCH -n 1000
#SBATCH -p normal
#SBATCH -t 00:30:00
#SBATCH --mail-user=dpederson@utexas.edu
#SBATCH --mail-type=end

mkdir strong_scaling
for p in 1 2 4 8 16 27 64 125 216 343 512 729 1000; do
	ibrun -np $p ../build/strong_oofdtd ;
	mkdir strong_scaling/np_$p ;
	mv proc*.log strong_scaling/np_$p/ ;
done
