#! /bin/bash

#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

/home/robin/bin/WHAM/WHAM.sh \
  -i Wallaby_heart_RRBS_rep1_Bock2023_GSM5853945_trimV1_bsmrk_macEug2_mrkdup.bam \
  -q 30 \
  -z /data/reference_genomes/macEug2/macEug2.sizes \
  -t 4 -l -m -C 4
