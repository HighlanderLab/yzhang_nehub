#!/bin/bash
#$ -N simulation
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
#$ -pe sharedmem 4
#$ -j y
#$ -V
#$ -cwd
#$ -o logs/sim_$TASK_ID.out

# Environment 
source conda
conda activate environment
module load bcftools
module load tabix
mkdir -p ../simulated4nehub

tasks=(
    "--name CON_SPLIT_MIG --path ../simulated4nehub/CON_SPLIT_MIG"
    "--name CON_SPLIT     --path ../simulated4nehub/CON_SPLIT"
    "--name CON           --path ../simulated4nehub/CON"
    "--name DEC_SPLIT_MIG --path ../simulated4nehub/DEC_SPLIT_MIG"
    "--name DEC_SPLIT     --path ../simulated4nehub/DEC_SPLIT"
    "--name DEC           --path ../simulated4nehub/DEC"
    "--length 8000000 --name CON --path ../simulated4nehub/CON_8Mb"
    "--length 8000000 --name DEC --path ../simulated4nehub/DEC_8Mb"
)
TASK=${tasks[$((SGE_TASK_ID - 1))]}
python sim.py $TASK

PATH_ARG=$(echo $TASK | grep -oP '(?<=--path )\S+')

for i in $(seq 0 9); do
    vcf="${PATH_ARG}/replicate${i}/sim_chr1.vcf"
    if [ -f "${vcf}" ]; then
        echo "  bgzip: ${vcf}"
        bgzip -f ${vcf}
        tabix -p -f vcf -f ${vcf}.gz
    elif [ -f "${vcf}.gz" ]; then
        echo "  already bgzipped: ${vcf}.gz"
        # 重新建 index
        tabix -p vcf -f ${vcf}.gz
    else
        echo "  WARNING: ${vcf} not found, skipping"
    fi
done

# generate pooled pop
base=path/to/save

for model in CON_SPLIT CON_SPLIT_MIG DEC_SPLIT DEC_SPLIT_MIG; do
  for r in replicate{0..9}; do
    d="$base/$model/$r"
    [ -f "$d/popLabel.txt" ] || { echo "skip $model/$r: no popLabel.txt"; continue; }
    awk 'BEGIN{FS="[ \t]+"; OFS="\t"}
         NR==1 { for (i=1;i<=NF;i++) if (tolower($i) ~ /^(population|group|pop)$/) col[i]=1
                 print; next }
         { for (i in col) $i="pop1"; print }' "$d/popLabel.txt" > "$d/popLabel_structure.txt"
    echo "$model/$r done"
  done
done