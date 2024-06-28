# JRA 2022
source /home/robin/sra2bw/sra2bw_functions.sh
CHROM_SIZES="/data/reference_genomes/mm10/mm10.sizes"
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"





for x in *CnT_NC*bam; do
	FILE=${x//.bam/}
	callPeaks narrow
done

