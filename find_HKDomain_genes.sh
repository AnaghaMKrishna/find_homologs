#! /usr/bin/bash
# usage: ./find_homologs.sh <query.faa> <subject.fna> <bedfile> <output file>

# Read args
HKDomain_query_file=$1
subject_file=$2
bed_file=$3
HKDomain_genes_file=$4

#Create temp files for storing intermediate outputs
#tblastn_output_temp=$(mktemp)
# tblastn command matches the query protein sequence with genome sequence in 
# subject files and outputs the standard 12 columns and query length.
# Use awk for column-based filtering to extract rows which match more than 
# 30% and have more than 90% match length and redirect to output file 
tblastn \
    -query $HKDomain_query_file \
    -subject $subject_file \
    -task tblastn \
    -outfmt '6 std qlen' \
| awk '{if ($3 > 30.000 && $4 > 0.9*$13 ) print $0;}' > tblastn_output_temp

#cut -f 1,9,10 tblastn_output_temp > tblastn_output_truncated_temp
awk 'FNR == NR {homologs_position[$2,$3]; next} {for (pos in homologs_position) {split(pos, r, SUBSEP); if ( r[1] >= $2 && r[2] <= $3) print $0}}' <(cut -f 1,9,10 tblastn_output_temp) $bed_file

#awk 'FNR==NR {a[$2]; next} { for (i in a) if (i > $2) print $0 }' tblastn_output_truncated_temp $bed_file
# awk 'FNR==NR {a[$2]; next} { for (i in a) if (i > $2) print $0}' tblastn_output_truncated_temp $bed_file
# awk 'FNR == NR {a[$2]; next} { for (i in a) print i}' tblastn_output_truncated_temp
# awk 'FNR == NR {a[$2]; b[$3]; next} { for (i in a) if (i > $2) print}' IFS='\t' tblastn_output_truncated_temp $bed_file
# awk 'FNR == NR {a[$2]; b[$3]; next} { for (i in a; j in b) if (i>$2 && j<$3) print}' tblastn_output_truncated_temp $bed_file
# awk 'FNR == NR {a[$2]; a[$3]; next} $3 in a' file1.txt file2.txt
#FNR==NR {a[$0];next} {for (i in a) if ($0~i) print}' file1 file2
# Count the number of lines in output file to obtain number of matches
#cat tblastn_output_temp 1163209