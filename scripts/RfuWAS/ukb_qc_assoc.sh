# UKB QC and association analysis
### 1. Sample QC
R qc.R

### 2. Variant QC
cd data/ukb/qc
plink2 --bfile ~/airr/data_from_xiaowei/ukb/bfiles_all_imputed/ukb_imp_chr6_v3.7.140998851.143510972 --set-all-var-ids @_#_\$r_\$a --new-id-max-allele-len 101 error --make-pgen --sort-vars --out ukb_imp_trb
plink2 --bfile ~/airr/data_from_xiaowei/ukb/bfiles_all_imputed/ukb_imp_chr6_v3.6.27477797.34448354 --set-all-var-ids @_#_\$r_\$a --new-id-max-allele-len 244 error --make-pgen --sort-vars --out ukb_imp_hla
plink2 --pfile ukb_imp_trb --pmerge ukb_imp_hla --make-pgen --out ukb_imp_trb_hla

plink2 --pfile ukb_imp_trb_hla --extract variants.pass.info.rsid --make-pgen --out ukb_imp_pass

plink2 --pfile ukb_imp_pass --geno 0.05 --make-pgen --out ukb_imp_missing

plink2 --pfile ukb_imp_missing --maf 0.001 --make-pgen --out ukb_imp_maf

plink2 --pfile ukb_imp_maf --hwe 1e-10 keep-fewhet --make-pgen --out ukb_imp_hwe

### 3. Phenotype data
### Phecodes
singularity exec \
   --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf \
   bioconductor_docker_latest.sif \
   bash -c 'cd /home2/s438989/airr/othercode/createUKBphenome && Rscript ./scripts/function.createUKBphenome.r'

R qc.R

##### END OF QC ##########################
##########################################

### 4. Fit lasso model
### 4.1 Find intersection between UKB and teddy
R lasso/map_ukb_teddy_variant.R

### 4.2 extract intersection variants from UKB and teddy
# Filename: trb.hla.ukb
cd ~/airr/data/ukb/lasso_variant_sets/
qc_dir=~/airr/data/ukb/qc/
snpset_label=".ukbflip"

## From UKB, extract variants in teddy, and rename to teddy name, and update allele to teddy allele. The ukb_sampleqc_variantqc_matrixeqtl.traw is for RFU prediction. We only use ID in RFU prediction, and ID has been converted to TEDDY name. ID, CHR and POS are in hg19. REF and ALT are in TEDDY order.
plink2 --pfile ${qc_dir}ukb_imp_hwe --extract ${qc_dir}ukb.varid --make-pgen --out ukb_intersect
plink2 --pfile ukb_intersect --update-name ${qc_dir}ukb.old.new.varid --make-pgen --out ukb_intersect_teddyname
plink2 --pfile ukb_intersect_teddyname --ref-allele ${qc_dir}ukb.ref 2 1 --make-pgen --out ukb_intersect_teddyname_teddyallele

plink2 --pfile ukb_intersect_teddyname_teddyallele --keep ${qc_dir}sample_filtered.id --export Av --out ukb_sampleqc_variantqc_Av
cut --complement -f1,3-6 ukb_sampleqc_variantqc_Av.traw > ukb_sampleqc_variantqc_matrixeqtl.traw

## From TEDDY, extract variants in UKB. The trb.hla.ukb file is for lasso model training on TEDDY dataset. ID is TEDDY name in hg19, CHR and POS are in hg38.
plink2 --pfile ../../gwas/sel_trb_hla_hg19 --extract ${qc_dir}teddy.varid --make-pgen --out trb.hla${snpset_label}
plink2 --pfile ../../gwas/sel_trb_hla_hg19 --extract ${qc_dir}teddy.varid --export Av --out trb.hla${snpset_label}_Av
cut --complement -f1,3-6 trb.hla${snpset_label}_Av.traw > trb.hla${snpset_label}_matrixeqtl.traw

plink2 --pfile ~/airr/data/ukb/lasso_variant_sets/trb.hla.ukbflip --freq --out ~/airr/data/ukb/lasso_variant_sets/trb.hla.ukbflip.freq

### 4.3 Extract variants for lasso model
declare -A chr_num from_bp to_bp suffixes
# For trb variants
chr_num[0]=7 from_bp[0]=141299011 to_bp[0]=143813287 suffixes[0]="trb${snpset_label}"
chr_num[1]=6 from_bp[1]=27510120 to_bp[1]=34480577 suffixes[1]="hla${snpset_label}"

# For trbonly
chr_num[2]=7 from_bp[2]=142299011 to_bp[2]=142813287 suffixes[2]="trbonly${snpset_label}"
chr_num[3]=6 from_bp[3]=28510120 to_bp[3]=33480577 suffixes[3]="hlaonly${snpset_label}"

# For trbsig
chr_num[4]=7 from_bp[4]=142243586 to_bp[4]=142931140 suffixes[4]="trbsig${snpset_label}"
chr_num[5]=6 from_bp[5]=30797802 to_bp[5]=33223718 suffixes[5]="hlasig${snpset_label}"

for ((i=0; i<${#chr_num[@]}; i++))
do
    prefix="${suffixes[$i]}"
    
    # Extract variant
    plink2 --pfile trb.hla${snpset_label} --chr ${chr_num[$i]} --from-bp ${from_bp[$i]} --to-bp ${to_bp[$i]} --export Av --out ${prefix}_Av
    cut --complement -f1,3-6 ${prefix}_Av.traw > ${prefix}_matrixeqtl.traw
done

Rscript ~/airr/scripts/lasso/rfu_lasso.R 1 trb${snpset_label} data/ukb/lasso_variant_sets/
Rscript ~/airr/scripts/lasso/rfu_lasso.R 1 trb${snpset_label}_hla${snpset_label} data/ukb/lasso_variant_sets/
Rscript ~/airr/scripts/lasso/rfu_lasso.R 1 trbonly${snpset_label} data/ukb/lasso_variant_sets/
Rscript ~/airr/scripts/lasso/rfu_lasso.R 1 trbonly${snpset_label}_hlaonly${snpset_label} data/ukb/lasso_variant_sets/
Rscript ~/airr/scripts/lasso/rfu_lasso.R 1 trbsig${snpset_label} data/ukb/lasso_variant_sets/
Rscript ~/airr/scripts/lasso/rfu_lasso.R 1 trbsig${snpset_label}_hlasig${snpset_label} data/ukb/lasso_variant_sets/

##### 4.4 GCTA-GREML
cd ~/airr/data/ukb/lasso_variant_sets/
mkdir gcta
# Set arrays
chr7=("141299011" "142299011" "142243586")
chr7_to=("143813287" "142813287" "142931140")

chr6=("27510120" "28510120" "30797802")
chr6_to=("34480577" "33480577" "33223718")

trbnames=("trb" "trbonly" "trbsig")
hlanames=("hla" "hlaonly" "hlasig")

# Loop over the arrays
for i in ${!trbnames[@]}; do
    plink2 --pfile trb.hla.ukbflip --chr 7 --from-bp ${chr7[$i]} --to-bp ${chr7_to[$i]} --merge-par --make-bed --out ukbflip_${trbnames[$i]}
    plink2 --pfile trb.hla.ukbflip --chr 6 --from-bp ${chr6[$i]} --to-bp ${chr6_to[$i]} --merge-par --make-bed --out ukbflip_${hlanames[$i]}
    plink2 --bfile ukbflip_${trbnames[$i]} --pmerge ukbflip_${hlanames[$i]}.bed ukbflip_${hlanames[$i]}.bim ukbflip_${hlanames[$i]}.fam --make-bed --out ukbflip_${trbnames[$i]}_${hlanames[$i]} --delete-pmerge-result --indiv-sort '0'
done

# For each command, execute gcta
for cmd in trb trb_hla trbonly trbonly_hlaonly trbsig trbsig_hlasig; do
    ~/airr/packages/gcta/gcta-1.94.1 --bfile ukbflip_$cmd --make-grm --out gcta/ukbflip_$cmd --thread-num 28
done

for cmd in hla hlaonly hlasig; do
    ~/airr/packages/gcta/gcta-1.94.1 --bfile ukbflip_$cmd --make-grm --out gcta/ukbflip_$cmd --thread-num 28
done

mkdir ~/airr/results/gcta/gcta_ukbf/
sbatch ~/airr/scripts/lasso/sbatch_gcta.sh

##### 4.5 Fit final model
R rfu_lasso_analyze.R
R fit_final_lasso.R

## 5. Predict RFU
for i in {1..12}; do
    sbatch ~/airr/scripts/lasso/predict_lasso_ukb.sh $i
done

time singularity exec \
--bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf \
bioconductor_docker_latest.sif Rscript /home2/s438989/airr/scripts/lasso/predict_lasso_ukb_merge.R &

## 6. Association
# The first parameter is phenotype file
# The second parameter is the prefix of the output file
sbatch sbatch_assoc.sh ~/airr/data/ukb/qc/UKB_PHENOME_QC.txt ukbflip_hg19

R assoc_ukb_rfu_analyze.R

sbatch sbatch_cancer.sh
R cancer_analyze.R
plink2 --pfile ~/airr/data/ukb/lasso_variant_sets/ukb_intersect_teddyname_teddyallele --keep ~/airr/data/ukb/stratify_tmp/phecode165.1_samples --extract  --export A --out ~/airr/data/ukb/stratify_tmp/phecode165.1_A

## 7. GSEA
mkdir ~/airr/results/interpretation/gsea_3827
gsea-cli.sh GSEA -res ~/airr/results/interpretation/gesa_expr.txt -cls ~/airr/results/interpretation/gesa_pheno_3827.cls -gmx ~/airr/packages/GSEA_Linux_4.3.2/h.all.v2023.2.Hs.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Signal2Noise -sort real -order descending -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out ~/airr/results/interpretation/gesa_3827

