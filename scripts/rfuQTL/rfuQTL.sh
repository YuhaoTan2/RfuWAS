# Codes for rfuQTL
##### Stats of bcf file
cd data
bcftools stats wgs_phased_variant_data.bcf > wgs_phased_variant_data.stats
plot-vcfstats wgs_phased_variant_data.stats -p vcfstats

##################################################
##### Begin of QC
##################################################
time bcftools norm -Ou -m - wgs_phased_variant_data.bcf | bcftools view -Ob -S wgs_sample_id.txt > selsamples.bcf # 30 min
time plink2 --bcf selsamples.bcf --double-id --allow-extra-chr --set-all-var-ids @_#_\$r_\$a --new-id-max-allele-len 101 error --split-par hg38 --make-pgen --update-sex selsamples_sex.tsv --out sel_plink2 # 4.5 min

##### 1. missing
# Delete SNPs with missingness >0.02, Delete individuals with missingness >0.02.
plink2 --pfile sel_plink2 --geno 0.02 --mind 0.02 --make-pgen --out sel_miss

##### 2. sex descrepancy
module add plink/1.9
plink2 --pfile sel_miss --merge-par --make-bed --chr X,Y --out sel_for_sex_check
plink --split-x hg38 --bfile sel_for_sex_check --make-bed --out sel_sex
plink --bfile sel_sex --check-sex --out sel_sex

grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
plink2 --bfile sel_miss --remove sex_discrepancy.txt --make-pgen --out sel_sex

##### 3. MAF
plink2 --pfile sel_miss --maf 0.05 --make-pgen --out sel_maf

##### 4. HWE
plink2 --pfile sel_maf --hwe 1e-6 keep-fewhet --make-pgen --out sel_hwe

##### 5. heterozygosity
plink2 --pfile sel_hwe --indep-pairwise 200 100 0.1 --out sel_indepSNP
plink2 --pfile sel_hwe --extract sel_indepSNP.prune.in --autosome --het --out sel_R_check

Rscript --no-save ../scripts/rfuQTL/heterozygosity_outliers_list.R

sed 's/"// g' sel_fail-het-qc.txt | awk '{print$1, $2}'> sel_het_fail_ind.txt

plink2 --pfile sel_hwe --remove sel_het_fail_ind.txt --make-pgen --out sel_hetero

##### 6. relatedness
plink2 --pfile sel_hetero --king-cutoff 0.0884 --make-pgen --out sel_rel
cat sel_rel.king.cutoff.out.id

##### 7. PCA
# LD pruning reduces the risk of getting PCs based on just a few genomic regions
plink2 --pfile sel_rel --extract sel_indepSNP.prune.in --autosome --pca --out sel_pca
plink2 --pfile sel_rel --extract sel_indepSNP.prune.in --autosome --export A-transpose --out sel_indepSNP_Av
# Predict ancestry
awk 'NR>1 {print $1}' sel_rel.psam > sel_rel_sample_id.txt
time bcftools view -Oz -S sel_rel_sample_id.txt wgs_phased_variant_data.bcf > sel_rel_whole_variants.vcf.gz 
~/airr/packages/GrafPop1.0/ExtractAncSnpsFromVcfGz.pl sel_rel_whole_variants.vcf.gz sel_rel_anc.vcf 
~/airr/packages/GrafPop1.0/grafpop sel_rel_anc.vcf sel_rel_pop.txt 
~/airr/packages/GrafPop1.0/PlotGrafPopResults.pl sel_rel_pop.txt sel_rel_pop.png
~/airr/packages/GrafPop1.0/SaveSamples.pl sel_rel_pop.txt sel_rel_smps.txt
# Check association between PC and ancestry
R gwas_qc.R

##################################################
##### SNPs QC ends, RFUs QC starts
##################################################

##### 8. TCR number > 2500
python3 gwas_pheno.py
R gwas_qc.R

##### 9. RFU normalize
python3 gwas_pheno.py

##### 10. PEER
R gwas_qc.R
nohup ~/airr/packages/install/bin/peertool -f gwas/subject_normalize_peer_tcrnum.tab -n 2 -o gwas/rfu_normalize_peer2_tcrnum --add_mean -c gwas/subject_rfu_peer_tcrnum.tab -i 1000  >> ../log/peer.out 2>&1 &


##### 11. MatrixEQTL Association
plink2 --pfile qctmp/sel_rel --keep tcrnum_sample.txt --make-pgen --out sel_tcrnum
plink2 --pfile sel_tcrnum --export Av --out sel_tcrnum_Av
cut --complement -f1,3-6 sel_tcrnum_Av.traw > sel_tcrnum_matrixeqtl.traw
sbatch ../scripts/rfuQTL/rfu_eqtl.sh data/gwas/sel_tcrnum_matrixeqtl.traw data/gwas/tcrnum_residuals.tsv 0 tcrnum_residuals 5e-8

R rfu_eqlt_analyze.R

##### 12. Convert to vcf and HLA imputation
cd ~/airr/data/reference_genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gzip -d hg38ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gzip -d hg38.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gzip -d hg19ToHg38.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz
gzip -d hg19.fa.gz

module add picard/2.10.3

time bcftools index qctmp/selsamples.bcf
time bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%ID\n' qctmp/selsamples.bcf --regions chr6:26000000-35000000 | awk '{gsub(/^chr/,"",$1); print}' > gwas/selsamples.old.new.tsv
plink2 --pfile gwas/sel_tcrnum --chr 6 --from-mb 26 --to-mb 35  --export vcf-4.2 id-delim=, --output-chr chr26 --out gwas/sel_tcrnum_vcf --update-name gwas/selsamples.old.new.tsv
CrossMap.py  vcf reference_genome/hg38ToHg19.over.chain  gwas/sel_tcrnum_vcf.vcf  reference_genome/hg19.fa  gwas/sel_tcrnum_liftover.vcf --no-comp-allele --chromid s
bcftools sort gwas/sel_tcrnum_liftover.vcf -Oz -o gwas/sel_tcrnum_liftover.vcf.gz

# Impute using https://imputationserver.sph.umich.edu/index.html#!run/imputationserver-hla
cd data/TEDDY/HLA_imputation
7z x -p"mI0RLvX8jXMoop" chr_6.zip
zcat chr6.dose.vcf.gz | sed 's/PASS;GENOTYPED/PASS/' | bcftools view  -i 'R2>.7' -Oz -o chr6.R20.7.dose.vcf.gz
plink2 --vcf chr6.R20.7.dose.vcf.gz dosage=HDS --id-delim , --make-pgen --out chr6.R20.7.dose # 48030 variants

# when you only test for HLA alleles and amino acids (modify as necessry)
cat chr6.R20.7.dose.pvar | grep -v "#" | grep -P 'HLA_[A-Za-z0-9]+\*\d{2}:\d{2,3}(?!\:\d{2})'  | cut -f3 >  test_markers.txt
cat chr6.R20.7.dose.pvar | grep -v "#" | grep -P 'HLA_[A-Za-z0-9]+\*\d{2}:\d{2,3}(?!\:\d{2})'  | awk '{print $3,"T"}' >  test_markers_alleles.txt # this is to make sure to output the dosages of the "presence" of the allele coded as T, but not the "absence" coded as A.

plink2 --vcf chr6.R20.7.dose.vcf.gz dosage=HDS --id-delim , --export Av --extract test_markers.txt --export-allele test_markers_alleles.txt --geno 1 --hwe 1e-4 --maf 0.01 --out chr6.R20.7.dose.sel # 59 variants
cut --complement -f1,3-6 chr6.R20.7.dose.sel.traw > chr6.R20.7.dose_matrixeqtl.traw
sbatch ~/airr/scripts/rfuQTL/rfu_eqtl.R data/TEDDY/HLA_imputation/chr6.R20.7.dose_matrixeqtl.traw data/gwas/tcrnum_residuals.tsv 0 hla_tcrnum_residuals 1
plink2 --vcf chr6.R20.7.dose.vcf.gz dosage=HDS --extract test_markers.txt --geno 1 --hwe 1e-4 --maf 0.01 --make-pgen --out chr6.R20.7.dose.sel
plink2 --pfile chr6.R20.7.dose.sel --freq --out chr6.R20.7.dose.sel

plink2 --vcf chr6.R20.7.dose.vcf.gz dosage=HDS --id-delim , --export Av --extract test_markers.txt --export-allele test_markers_alleles.txt --keep ~/airr/data/gwas/tcr1e4_sample.txt --geno 1 --hwe 1e-4 --maf 0.01 --out chr6.R20.7.dose.tcr1e4 # 61 variants, 318 samples
cut --complement -f1,3-6 chr6.R20.7.dose.tcr1e4.traw > chr6.R20.7.dose.tcr1e4_matrixeqtl.traw

