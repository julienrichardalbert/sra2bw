#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

for x in /scratch/WHAM/bams/Mouse*bam; do

  /home/robin/bin/WHAM/WHAM.sh \
    -i $x \
    -q 30 \
    -z /data/reference_genomes/mm10/mm10.sizes \
    -t 4 -m -C 4

done
