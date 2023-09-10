#! /usr/bin/bash
# usage: ./find_homologs.sh <query file> <subject file> <output file>

# tblastn command matches the query protein sequence with genome sequence in 
# subject files and outputs the standard 12 columns and query length.
# Use awk for column-based filtering to extract rows which match more than 
# 30% and have more than 90% match length and redirect to output file 
tblastn -query $1 -subject $2 -task tblastn -outfmt '6 std qlen' | 
awk '{if ($3 > 30.000 && $4 > 0.9*$13 ) print $0;}' > $3

# Count the number of lines in output file to obtain number of matches
wc -l $3 