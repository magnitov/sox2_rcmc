DATA_PATH=/DATA/projects/bvs_alleles/sox2_rcmc_custom

mkdir -p ${DATA_PATH}/pairs
mkdir -p ${DATA_PATH}/stats

for SAMPLE in CM1846 CM2167 CM2267 CM2094 CM2266 CM2291
do
	if [ "$SAMPLE" == "CM1846" ] || [ "$SAMPLE" == "CM2094" ]; then
		for REPLICATE in rep1 rep2_batch1 rep2_batch2
		do
			pairtools parse --min-mapq 0 --add-columns mapq,XA,NM,AS,XS --drop-sam --walks-policy mask \
				-c ${DATA_PATH}/custom_genomes/CM_129_CAST_custom.chrom.sizes \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.unphased.XA.pairs.gz \
				${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}.mapped.XA.bam
            
			pairtools phase --phase-suffixes _129 _CAST --tag-mode XA \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.XA.pairs.gz \
				${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.unphased.XA.pairs.gz
            
			pairtools sort --nproc 32 \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.XA.pairs.gz \
				${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.XA.pairs.gz
		done

		pairtools dedup --mark-dups --backend cython --max-mismatch 1 \
			--extra-col-pair "phase1" "phase2" --output-dups - --output-unmapped - \
			--output-stats ${DATA_PATH}/stats/${SAMPLE}_rep1.phased.sorted.merged.XA.dedup.stats \
			-o ${DATA_PATH}/pairs/${SAMPLE}_rep1.phased.sorted.merged.XA.nodup.pairs.gz \
			${DATA_PATH}/pairs/${SAMPLE}_rep1.phased.sorted.XA.pairs.gz
                
		pairtools merge --nproc 32 \
			-o ${DATA_PATH}/pairs/${SAMPLE}_rep2.phased.sorted.merged.XA.pairs.gz \
			${DATA_PATH}/pairs/${SAMPLE}_rep2_batch*.phased.sorted.XA.pairs.gz
        
		pairtools dedup --mark-dups --backend cython --max-mismatch 1 \
			--extra-col-pair "phase1" "phase2" --output-dups - --output-unmapped - \
			--output-stats ${DATA_PATH}/stats/${SAMPLE}_rep2.phased.sorted.merged.XA.dedup.stats \
			-o ${DATA_PATH}/pairs/${SAMPLE}_rep2.phased.sorted.merged.XA.nodup.pairs.gz \
			${DATA_PATH}/pairs/${SAMPLE}_rep2.phased.sorted.merged.XA.pairs.gz
	fi
    
	if [ "$SAMPLE" == "CM2167" ] || [ "$SAMPLE" == "CM2267" ] || [ "$SAMPLE" == "CM2266" ] || [ "$SAMPLE" == "CM2291" ]; then
		for REPLICATE in rep1_batch1 rep1_batch2 rep2_batch1 rep2_batch2
		do
			pairtools parse --min-mapq 0 --add-columns mapq,XA,NM,AS,XS --drop-sam --walks-policy mask \
				-c ${DATA_PATH}/custom_genomes/CM_129_CAST_custom.chrom.sizes \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.unphased.XA.pairs.gz \
				${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}.mapped.XA.bam
            
			pairtools phase --phase-suffixes _129 _CAST --tag-mode XA \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.XA.pairs.gz \
				${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.unphased.XA.pairs.gz
            
			pairtools sort --nproc 32 \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.XA.pairs.gz \
				${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.XA.pairs.gz
		done

		for REPLICATE in rep1 rep2
		do
			pairtools merge --nproc 32 \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.merged.XA.pairs.gz \
				${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}_batch*.phased.sorted.XA.pairs.gz
        
			pairtools dedup --mark-dups --backend cython --max-mismatch 1 \
				--extra-col-pair "phase1" "phase2" --output-dups - --output-unmapped - \
				--output-stats ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.phased.sorted.merged.XA.dedup.stats \
				-o ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.merged.XA.nodup.pairs.gz \
				${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.merged.XA.pairs.gz
		done
	fi
done
