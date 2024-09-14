#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/early_SC2_trajectory/gisaid_scripts
#$ -S /bin/bash
#$ -pe smp 24
#$ -j y
#$ -R y
#$ -N generate

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla /SAN/ugi/HAP_VAP/early_SC2_trajectory/gisaid_scripts/generate_presence_matrix.cs.R
