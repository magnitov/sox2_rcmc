DATA_PATH=/DATA/projects/bvs_alleles/sox2_rcmc_custom

mkdir -p ${DATA_PATH}/coolers

for SAMPLE in CM1846 CM2167 CM2267 CM2094 CM2266 CM2291
do
	for REPLICATE in rep1 rep2
	do
		cooler cload pairs --assembly mm10_129S1_CAST -c1 2 -p1 3 -c2 4 -p2 5 \
			${DATA_PATH}/custom_genomes/CM_129_CAST_custom_merged.chrom.sizes:1000 \
			${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs \
			${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.1000.cool

		cooler zoomify --resolutions 1000,2000,5000,10000,20000,50000,100000,200000,500000 \
			--balance -o ${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.mcool \
			${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.1000.cool
            
		for RULE in ds ss
		do
			cooler cload pairs --assembly mm10_129S1_CAST -c1 2 -p1 3 -c2 4 -p2 5 \
				${DATA_PATH}/custom_genomes/CM_129_CAST_custom_merged.chrom.sizes:1000 \
				${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.pairs \
				${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.1000.cool

			cooler zoomify --resolutions 1000,2000,5000,10000,20000,50000,100000,200000,500000 \
				--balance -o ${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.mcool \
				${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.1000.cool
            
			cooler cload pairs --assembly mm10_129S1_CAST -c1 2 -p1 3 -c2 4 -p2 5 \
				${DATA_PATH}/custom_genomes/CM_129_CAST_custom_merged.chrom.sizes:1000 \
				${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.pairs \
				${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.1000.cool

			cooler zoomify --resolutions 1000,2000,5000,10000,20000,50000,100000,200000,500000 \
				--balance -o ${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.mcool \
				${DATA_PATH}/coolers/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.1000.cool
		done
	done
    
	cooler merge ${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.1000.cool \
		${DATA_PATH}/coolers/${SAMPLE}_rep1.sox2_rcmc.1000.cool \
		${DATA_PATH}/coolers/${SAMPLE}_rep2.sox2_rcmc.1000.cool
    
	cooler zoomify --resolutions 1000,2000,5000,10000,20000,50000,100000,200000,500000 \
		--balance -o ${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.mcool \
		${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.1000.cool

	for RULE in ds ss
	do
		cooler merge ${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.129_allele.${RULE}.1000.cool \
			${DATA_PATH}/coolers/${SAMPLE}_rep1.sox2_rcmc.129_allele.${RULE}.1000.cool \
			${DATA_PATH}/coolers/${SAMPLE}_rep2.sox2_rcmc.129_allele.${RULE}.1000.cool
    
		cooler zoomify --resolutions 1000,2000,5000,10000,20000,50000,100000,200000,500000 \
			--balance -o ${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.129_allele.${RULE}.mcool \
			${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.129_allele.${RULE}.1000.cool
        
		cooler merge ${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.CAST_allele.${RULE}.1000.cool \
			${DATA_PATH}/coolers/${SAMPLE}_rep1.sox2_rcmc.CAST_allele.${RULE}.1000.cool \
			${DATA_PATH}/coolers/${SAMPLE}_rep2.sox2_rcmc.CAST_allele.${RULE}.1000.cool
    
		cooler zoomify --resolutions 1000,2000,5000,10000,20000,50000,100000,200000,500000 \
			--balance -o ${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.CAST_allele.${RULE}.mcool \
			${DATA_PATH}/coolers/${SAMPLE}.sox2_rcmc.CAST_allele.${RULE}.1000.cool
	done

	rm ${DATA_PATH}/coolers/${SAMPLE}*.cool
done
