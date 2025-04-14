### Download VCF with SNPs
wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz
wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz

### Filter only homozygous SNPs
java -jar /DATA/users/m.magnitov/software/snpEff_v4_3p/SnpSift.jar filter -v "(isHom( GEN[0] ))" 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz > 129S1_SvImJ.mgp.v5.snps.dbSNP142.homozygous.vcf
java -jar /DATA/users/m.magnitov/software/snpEff_v4_3p/SnpSift.jar filter -v "(isHom( GEN[0] ))" CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz > CAST_EiJ.mgp.v5.snps.dbSNP142.homozygous.vcf
bgzip 129S1_SvImJ.mgp.v5.snps.dbSNP142.homozygous.vcf
bgzip CAST_EiJ.mgp.v5.snps.dbSNP142.homozygous.vcf
bcftools index 129S1_SvImJ.mgp.v5.snps.dbSNP142.homozygous.vcf.gz
bcftools index CAST_EiJ.mgp.v5.snps.dbSNP142.homozygous.vcf.gz

### Download genome
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna//Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
zcat Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz | bgzip -c > Mus_musculus.GRCm38.dna_sm.toplevel.fa.bgz
samtools faidx Mus_musculus.GRCm38.dna_sm.toplevel.fa.bgz
rm Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz

### Incorporate SNPs into the genome
bcftools consensus --fasta-ref Mus_musculus.GRCm38.dna_sm.toplevel.fa.bgz \
	--haplotype 1 129S1_SvImJ.mgp.v5.snps.dbSNP142.homozygous.vcf.gz | sed -E 's/(>[^[:space:]]+).*/\1_129/g' | bgzip -c > GRCm38_129_snpsonly.fa.gz
bcftools consensus --fasta-ref Mus_musculus.GRCm38.dna_sm.toplevel.fa.bgz \
	--haplotype 1 CAST_EiJ.mgp.v5.snps.dbSNP142.homozygous.vcf.gz | sed -E 's/(>[^[:space:]]+).*/\1_CAST/g' | bgzip -c > GRCm38_CAST_snpsonly.fa.gz

# Merge and index combined genome
zcat GRCm38_129_snpsonly.fa.gz GRCm38_CAST_snpsonly.fa.gz > GRCm38_129_CAST_snpsonly.fa
bwa index GRCm38_129_CAST_snpsonly.fa
