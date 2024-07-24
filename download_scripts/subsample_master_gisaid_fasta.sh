large=../data/genomes/gisaid/gisaid_sequences.230324.renamed.fna
id_list=../data/metadata/gisaid_metatadata.220324.1Apr20.accessions_only.txt
out=../data/genomes/gisaid/gisaid_sequences.230324.1Apr20.fna

seqkit split2 -j 20 -p 100 $large -O tmp


#seqtk subseq \
#    $large \
#    $id_list \
#    > $out
