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

cut -f 1,9,10 $tblastn_output_temp

# Count the number of lines in output file to obtain number of matches
cat $tblastn_output_temp