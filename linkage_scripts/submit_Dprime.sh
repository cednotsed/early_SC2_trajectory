#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/early_SC2_trajectory/linkage_scripts
#$ -S /bin/bash
#$ -pe smp 16
#$ -j y
#$ -R y
#$ -N Dprime

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate r_env

dataset=pirola

Rscript --vanilla /SAN/ugi/HAP_VAP/early_SC2_trajectory/linkage_scripts/calculate_Dprime.cs.R $dataset
