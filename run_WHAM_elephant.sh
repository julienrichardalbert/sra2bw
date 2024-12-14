#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

/home/robin/bin/WHAM/WHAM.sh \
  -i Elephant_heart_RRBS_rep1_Bock2023_GSM5853042_trimV1_bsmrk_loxAfr3_mrkdup.bam \
  -q 30 \
  -z /data/reference_genomes/loxAfr3/loxAfr3.sizes \
  -t 4 -l -m
