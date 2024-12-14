#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

/home/robin/bin/WHAM/WHAM.sh \
  -i /data/bams/Armadillo_heart_RRBS_rep1_Bock2023_GSM5854456_trimV1_bsmrk_dasNov3_mrkdup.bam \
  -b dasNov3_CGIs.bed \
  -q 30 \
  -z dasNov3.sizes \
  -t 4 -l -m -C 4
