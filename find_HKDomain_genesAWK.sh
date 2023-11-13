#! /usr/bin/bash
# usage: ./find_homologs.sh <query.faa> <subject.fna> <bedfile> <output file>

# Read args
HKDomain_query_file=$1
subject_file=$2
bed_file=$3
HKDomain_genes_file=$4

#Create temp files for storing intermediate outputs
tblastn_output_temp=$(mktemp)

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

#head $tblastn_output_temp

# Loop over start and end positions of homologs and compare if it lies with 
# the boundary of any genes in BED file and filter for unique gene names
awk 'FNR == NR {
        homologs_position[$2,$3]; 
        next
    } 
    {
        for (pos in homologs_position) 
        {
            split(pos, pos_arr, SUBSEP); 
            if ( pos_arr[1] > $2 && pos_arr[2] <= $3)
                print $4
        }
    }
    ' <(cut -f 1,9,10 $tblastn_output_temp) $bed_file | sort | uniq > $HKDomain_genes_file

# Get the organism name from input file name
org=${subject_file##*/}
echo "The number of homologous HK domain genes identified for ${org%.*} is $(wc -l $HKDomain_genes_file | cut -d' ' -f1) and written to file $(wc -l $HKDomain_genes_file | cut -d' ' -f2)"

# Clean up temp file
rm $tblastn_output_temp