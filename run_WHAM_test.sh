#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

/home/robin/bin/WHAM/WHAM.sh \
  -i Mouse_heart_RRBS_rep1_Bock2023_GSM5852816_trimV1_bsmrk_mm10_mrkdup.bam \
  -q 30 \
  -z /data/reference_genomes/mm10/mm10.sizes \
  -t 4 -m -C 4
