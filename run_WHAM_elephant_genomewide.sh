#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

for x in /scratch/WHAM/bams/Elep*bam; do

  /home/robin/bin/WHAM/WHAM.sh \
    -i $x \
    -q 30 \
    -z /data/reference_genomes/xenarthra/elephant/data/GCF_030014295.1/GCF_030014295.1_mLoxAfr1.hap2_genomic.sizes \
    -t 4 -m -C 4

done


