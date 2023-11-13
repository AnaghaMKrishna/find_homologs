#! /usr/bin/bash
# usage: ./find_homologs.sh <query.faa> <subject.fna> <bedfile> <output file>

# Read args
HKDomain_query_file=$1
subject_file=$2
bed_file=$3
HKDomain_genes_file=$4

#Create temp files for storing intermediate outputs
tblastn_output_temp=$(mktemp)
tblastn_output_trimmed=$(mktemp)
gene_names_temp=$(mktemp)

# tblastn command matches the query protein sequence with genome sequence in 
# subject files and outputs the standard 12 columns and query length.
# Use awk for column-based filtering to extract rows which match more than 
# 30% and have more than 90% match length and redirect to output file 
tblastn \
    -query $HKDomain_query_file \
    -subject $subject_file \
    -task tblastn \
    -outfmt '6 std qlen' \
| awk '{if ($3 > 30.000 && $4 > 0.9*$13 ) print $0;}' > $tblastn_output_temp

cut -f 1,9,10 $tblastn_output_temp > $tblastn_output_trimmed
while read id sstart send 
do
    while read bed_id bed_start bed_stop gene dot direction
    do
        if [[ $bed_start -lt $sstart && $send -le $bed_stop ]] ; then 
            echo $gene >> $gene_names_temp
        fi
    done < $bed_file
done < $tblastn_output_trimmed

# Get only unique gene names
sort $gene_names_temp | uniq > $HKDomain_genes_file

# Get the organism name from input file name
org=${subject_file##*/}
echo "The number of homologous HK domain genes identified for ${org%.*} is $(wc -l $HKDomain_genes_file | cut -d' ' -f1) and written to file $(wc -l $HKDomain_genes_file | cut -d' ' -f2)"

# Clean up temp file
rm $tblastn_output_temp
rm $tblastn_output_trimmed
rm $gene_names_temp

# for orgs in $(ls -D ./);  do ../find_homologs/find_HKDomain_genes.sh ../HK_domain_Q2.faa 
# $(ls $orgs/*.fna) $(ls $orgs/*.bed) $(ls $orgs/*.bed | cut -d. -f1).txt; done