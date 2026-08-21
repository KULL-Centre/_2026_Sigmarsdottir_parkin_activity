#!/bin/bash

function join_by { local IFS="$1"; shift; echo "$*"; }

[[ -f samples.csv ]] && { echo "File samples.csv exists - cowardly refusing to overwrite!"; exit 2; } 

lib_name="prkn_bc"
lib_file="barcode_map.csv"

echo "# Input for match_and_merge.r script by kristoffer.johansson@bio.ku.dk" >> samples.csv
echo "# The output will be a data frame per library with a column per row in the samples section; i.e. per selection-replicate pair containting read counts per library member summed for technical replicates"  >> samples.csv
echo "" >> samples.csv

echo "SETTINGS" >> samples.csv
echo "# Section to give pairs of options and values" >> samples.csv
# echo "library_is_dna,true" >> samples.csv
# echo "library_is_protein,true" >> samples.csv
echo "library_is_barcode,true" >> samples.csv
echo "" >> samples.csv

echo "LIBRARIES" >> samples.csv
echo "# Section to list libraries that samples should be matched to (if any). Include control sequences using the name 'controls'" >> samples.csv
echo "# Library name,Library path" >> samples.csv
echo "$lib_name,$lib_file" >> samples.csv
echo "" >> samples.csv

echo "SAMPLES" >> samples.csv
echo "# Section to list sample files. All technical replicates should be given in same row" >> samples.csv
echo "# Library name,Selection,Replicate,Description,Raw file 1,Raw file 2,Raw file 3,Raw file 4,Raw file 5,Raw file 6,Raw file 7,Raw file 8,More if needed" >> samples.csv

# 555_S7_L003_counts_all.txt.gz
rep=1
sel="ao"
for bin in 1 2 3 4; do
    fields="$lib_name ${sel}-bin${bin} $rep none $(ls counts/55${bin}_S*_counts_all.txt.gz)"
    join_by , $fields >> samples.csv
done
sel="dmso"
for bin in 5 6 7 8; do
    fields="$lib_name ${sel}-bin$((bin-4)) $rep none $(ls counts/55${bin}_S*_counts_all.txt.gz)"
    join_by , $fields >> samples.csv
done

# 566-I_S16_L002_counts_all.txt.gz
rep=2
sel="ao"
for bin in 1 2 3 4; do
    fields="$lib_name ${sel}-bin${bin} $rep none $(ls counts/56${bin}-I_S*_counts_all.txt.gz)"
    join_by , $fields >> samples.csv
done
sel="dmso"
for bin in 5 6 7 8; do
    fields="$lib_name ${sel}-bin$((bin-4)) $rep none $(ls counts/56${bin}-I_S*_counts_all.txt.gz)"
    join_by , $fields >> samples.csv
done

# 561-II_S19_L004_counts_all.txt.gz
rep=3
sel="ao"
for bin in 1 2 3 4; do
    fields="$lib_name ${sel}-bin${bin} $rep none $(ls counts/56${bin}-II_S*_counts_all.txt.gz)"
    join_by , $fields >> samples.csv
done
sel="dmso"
for bin in 5 6 7 8; do
    fields="$lib_name ${sel}-bin$((bin-4)) $rep none $(ls counts/56${bin}-II_S*_counts_all.txt.gz)"
    join_by , $fields >> samples.csv
done

