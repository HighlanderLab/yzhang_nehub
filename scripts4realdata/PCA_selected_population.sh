#!/bin/bash
#$ -N pca_sne
#$ -cwd
#$ -R y
#$ -pe sharedmem 8
#$ -l h_vmem=16G
#$ -l h_rt=05:25:00
#$ -j y
#$ -V
#$ -m ea

# Environment 
module load plink/2.0-dev     
# Load GCTA
module load gcta/1.94.1
input_PREFIX="path/to/vcf/file"
input_PostFIX=".vcf.gz"
output_PREFIX="path/to/save"
merge_list="path/to/merge/list"


# Clean merge-list.txt
> $merge_list

# Convert VCFs to PLINK format
for ((i = 1; i <= 18; i++)); do
    filename="${input_PREFIX}${i}${input_PostFIX}"
    output="${output_PREFIX}${i}"
    plink2 --vcf $filename --make-pgen --out $output
    echo $output >> $merge_list
done

# Merge all binary files
plink2 --pmerge-list $merge_list --make-bed --out ../data/merged


# Make GRM
gcta64 --bfile ../data02/merged --make-grm --out ../data/merged

# PCA for all individuals
gcta64 --grm ../data02/merged --pca --out ../data/merged_total

# PCA for subsets
sed 's/^/0\t/' ../data02/total.txt > ../data/total_new.txt
sed 's/^/0\t/' ../data02/Yorkshire.txt > ../data/Yorkshire_new.txt
sed 's/^/0\t/' ../data02/Landrace.txt > ../data/Landrace_new.txt
sed 's/^/0\t/' ../data02/Duroc.txt > ../data/Duroc_new.txt

gcta64 --grm ../data02/merged --pca --keep ../data/total_new.txt     --thread-num 6 --out ../data02/merged_total
gcta64 --grm ../data02/merged --pca --keep ../data/Yorkshire_new.txt --thread-num 6 --out ../data02/merged_Yorkshire
gcta64 --grm ../data02/merged --pca --keep ../data/Landrace_new.txt  --thread-num 6 --out ../data02/merged_Landrace
gcta64 --grm ../data02/merged --pca --keep ../data/Duroc_new.txt     --thread-num 6 --out ../data02/merged_Duroc






