DATA_PATH=/DATA/projects/bvs_alleles/sox2_rcmc_custom

mkdir -p ${DATA_PATH}/pairs_rcmc

for SAMPLE in CM1846 CM2167 CM2267 CM2094 CM2266 CM2291
do
	for REPLICATE in rep1 rep2
	do
		# Extract reads from Sox2 capture locus
		zcat ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.merged.XA.nodup.pairs.gz | grep '#' \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header
		zcat ${DATA_PATH}/pairs/${SAMPLE}_${REPLICATE}.phased.sorted.merged.XA.nodup.pairs.gz | grep -v '#' |\
			awk -F'\t' '{ if (($2=="3") && ($4=="3") && ($3>=33750000) && ($3<=35659386) && ($5>=33750000) && ($5<=35659386) && ($8!="DD")) print $0 }' \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.body
		cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header \
			${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.body \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs
		pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs \
			> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.stats
            
		rm ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.header ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.body

		#  Split reads into 129 and CAST alleles
		# DS: reads with both sides phased
		#
		grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header
		grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs |\
			awk -F'\t' '{ if ($19=="0" && $20=="0") print $0 }' \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.body
		grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs |\
			awk -F'\t' '{ if ($19=="1" && $20=="1") print $0 }' \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.body
                
		cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header \
			${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.body \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.ds.pairs
		pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.ds.pairs \
			> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.ds.stats
		cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header \
			${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.body \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.ds.pairs
		pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.ds.pairs \
			> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.ds.stats
                
		rm ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.header ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.body
        
		# SS: reads with at least on side phased
		#
		grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header
		grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs |\
			awk -F'\t' '{ if (($19=="0" && $20=="0") || ($19=="\." && $20=="0") || ($19=="0" && $20=="\.")) print $0 }' \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.body
		grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.pairs |\
			awk -F'\t' '{ if (($19=="1" && $20=="1") || ($19=="\." && $20=="1") || ($19=="1" && $20=="\.")) print $0 }' \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.body
                
		cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header \
			${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.body \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.ss.pairs
		pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.ss.pairs \
			> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.ss.stats
		cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.header \
			${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.body \
			> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.ss.pairs
		pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.ss.pairs \
			> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.ss.stats
                
		rm ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.header ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.body
	done
done


for SAMPLE in CM1846 CM2167 CM2267 CM2094 CM2266 CM2291
do
	for REPLICATE in rep1 rep2
	do
		for RULE in ds ss
		do      
			if [[ "$SAMPLE" == "CM1846" ]]; then
				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body

				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body
			fi
            

			if [[ "$SAMPLE" == "CM2167" ]]; then
				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body

				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body
			fi
            
            
			if [[ "$SAMPLE" == "CM2267" ]]; then
				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34656323) || ($3>34661999)) && (($5<34656323) || ($5>34661999))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body

				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body
			fi

            
			if [[ "$SAMPLE" == "CM2094" ]]; then
				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34601113) || ($3>34602127)) && (($5<34601113) || ($5>34602127))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body

				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs |\
					awk -F'\t' '{ if (($13==0) || ($14==0)) print $0 }' |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body
			fi
            
            
			if [[ "$SAMPLE" == "CM2266" ]]; then
				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34601113) || ($3>34602127)) && (($5<34601113) || ($5>34602127))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34656323) || ($3>34661999)) && (($5<34656323) || ($5>34661999))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body

				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body
			fi
            

			if [[ "$SAMPLE" == "CM2291" ]]; then
				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34656323) || ($3>34661999)) && (($5<34656323) || ($5>34661999))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body

				grep '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header
				grep -v '#' ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.pairs |\
					awk -F'\t' '{ if ((($3<34598480) || ($3>34603292)) && (($5<34598480) || ($5>34603292))) print $0 }' |\
					awk -F'\t' '{ if ((($3<34648847) || ($3>34652698)) && (($5<34648847) || ($5>34652698))) print $0 }' \
					> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body
			fi
        
			cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.header \
				${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.body \
				> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.pairs
			pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.pairs \
				> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.129_allele.${RULE}.filtered.stats
			cat ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.header \
				${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.body \
				> ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.pairs
			pairtools stats ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.pairs \
				> ${DATA_PATH}/stats/${SAMPLE}_${REPLICATE}.sox2_rcmc.CAST_allele.${RULE}.filtered.stats
		done
        
		rm ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.header ${DATA_PATH}/pairs_rcmc/${SAMPLE}_${REPLICATE}*.body
	done
done
