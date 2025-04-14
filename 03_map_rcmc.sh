DATA_PATH=/DATA/projects/bvs_alleles/sox2_rcmc_custom/

mkdir -p ${DATA_PATH}_custom/bam


for SAMPLE in CM1846 CM2094
do
	for REPLICATE in rep1 rep2_batch1 rep2_batch2
	do
		bwa mem -SP -t 24 ${DATA_PATH}/custom_genomes/${SAMPLE}_129_CAST_custom.fa \
			${DATA_PATH}/fastq/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${SAMPLE}_${REPLICATE}_R2.fastq.gz |\
			samtools view -@ 24 -b > ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}.mapped.XA.bam
	done
done

for SAMPLE in CM2167 CM2267 CM2266 CM2291
do
	for REPLICATE in rep1_batch1 rep1_batch2 rep2_batch1 rep2_batch2
	do
		bwa mem -SP -t 24 ${DATA_PATH}/custom_genomes/${SAMPLE}_129_CAST_custom.fa \
			${DATA_PATH}/fastq/${SAMPLE}_${REPLICATE}_R1.fastq.gz ${DATA_PATH}/fastq/${SAMPLE}_${REPLICATE}_R2.fastq.gz |\
			samtools view -@ 24 -b > ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}.mapped.XA.bam
	done
done
