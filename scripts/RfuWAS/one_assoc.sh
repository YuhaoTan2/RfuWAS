#!/bin/bash
phenotype_file=$2
# Read the first line of a file into a variable
exec 3<$2
read -r first_line <&3
exec 3<&-

# Convert tab-separated list to array
IFS=$'\t' read -ra elements <<< "$first_line"

# Get the total number of elements in the array and exclude the last element
total_elements=$((${#elements[@]} - 1))

# Calculate the number of elements per part
elements_per_part=$((total_elements / 14 + 1))

# Calculate the starting and ending index for the current part
start_index=$((($1 - 1) * elements_per_part))
end_index=$(($1 * elements_per_part - 1))

# Ensure end_index doesn't exceed the total number of elements
end_index=$(( end_index < total_elements-1 ? end_index : total_elements-1 ))

echo "start" $start_index "end" $end_index

# Process the sub-list for the current part
for ((i=start_index; i<=end_index; i++)); do
    element=${elements[$i]}
    echo $element

    python3 ~/airr/othercode/MetaXcan/software/PrediXcanAssociation.py \
    --expression_file ~/airr/results/lasso_ukbf/bestmodel.ukbflip_ukb_hg19_merge.txt \
    --input_phenos_file ~/airr/data/ukb/phenotype/${element}_permuted.tsv \
    --input_phenos_column ${element} \
    --covariates_file ~/airr/data/ukb/qc/sample_qc.kt \
    --covariates isMale PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 age \
    --output ~/airr/results/assoc_pheno/${3}/assoc_${element}_permuted.tsv \
    --mode logistic
done

