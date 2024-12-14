#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

for x in /scratch/WHAM/bams/Sloth*bam; do

  /home/robin/bin/WHAM/WHAM.sh \
    -i $x \
    -q 30 \
    -z /data/reference_genomes/xenarthra/sloth/data/GCF_015220235.1/GCF_015220235.1_mChoDid1.pri_genomic.sizes \
    -t 4 -m -C 4

done


