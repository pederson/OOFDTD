#!/bin/bash

#SBATCH -J oofdtd_weak
#SBATCH -o oofdtd_weak.o%j
#SBATCH -n 1000
#SBATCH -p normal
#SBATCH -t 00:30:00
#SBATCH --mail-user=dpederson@utexas.edu
#SBATCH --mail-type=end

mkdir weak_scaling
for p in 1 2 4 8 16 27 64 125 216 343 512 729 1000; do
	ibrun -np $p ../build/weak_oofdtd 100000 ;
	mkdir weak_scaling/np_$p ;
	mv proc*.log weak_scaling/np_$p/ ;
done
