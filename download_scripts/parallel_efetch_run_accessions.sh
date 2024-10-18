#acc_list=../data/metadata/sra_metadata/filtered_sra_accessions.txt
#out_path=/mnt/c/git_repos/wuhu_rooting/data/metadata/sra_metadata/runinfo.csv

#acc_list=../data/metadata/sra_metadata/filtered_sra_accessions.missing.txt
#out_path=/mnt/c/git_repos/wuhu_rooting/data/metadata/sra_metadata/runinfo.missing.csv

export acc_list=../data/metadata/sra_metadata/filtered_sra_accessions.delta.txt
export out_dir=../data/metadata/sra_metadata/runinfo_temp.delta

export NCBI_API_KEY=

cat $acc_list|xargs -P10 -I{} sh -c \
    'acc={}; \
        esearch \
            -db sra \
            -query $acc \
            | efetch \
            -format runinfo \
            |tail -n +2 \
            > $out_dir/$acc.csv'

