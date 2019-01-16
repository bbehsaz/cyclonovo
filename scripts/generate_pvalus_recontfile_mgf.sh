allreconstfile=$1

allspectfile=$2

output=$3

cyclonovoDir=$4

err=$5
threshcore=11

cat $allreconstfile | awk '($(NF)>11){print $(NF-3)}' | sort | uniq > "$output"_masses.txt

while read line; do grep "$line" "$allreconstfile" | awk '(NF>3 && $(NF-2)>11){n++; if (n<100){print $1}}' | sed 's/,/ /g' > "$output"_"$line"_reconst.txt  ; done < "$output"_masses.txt
while read line; do grep "$line" "$allreconstfile" | awk '(NF>3 && $(NF-2)>11){n++; if (n<100){print $0}}' > "$output"_"$line"_allinfo.txt  ; done < "$output"_masses.txt

cat $allspectfile  | python "$cyclonovoDir"/scripts/get_by_pepmass.py "$output"_masses.txt "$output"_


while read line; do if [ -f "$output"_"$line".mgf ]; then print_score "$output"_"$line".mgf "$output"_"$line"_reconst.txt --mass_seq_in --no_filter --product_ion_thresh "$err" --no_merge --make_cyclic --multiple_seq_file --concise_output | grep -v "^#" > "$output"_"$line"_pvals.txt; fi; done < "$output"_masses.txt 

while read line; do if [ -f "$output"_"$line".mgf ]; then paste "$output"_"$line"_allinfo.txt "$output"_"$line"_pvals.txt | awk '{print $(NF-1),$0}' | sort -nr | cut -f2- -d' ' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7}'  > "$output"_"$line"_reconst_pvals.txt; fi; done < "$output"_masses.txt 

while read line; do if [ -f "$output"_"$line"_allinfo.txt ]; then rm "$output"_"$line"_allinfo.txt; fi; done < "$output"_masses.txt

while read line; do if [ -f "$output"_"$line"_pvals.txt ]; then rm "$output"_"$line"_pvals.txt; fi; done < "$output"_masses.txt

while read line; do if [ -f  "$output"_"$line".mgf ]; then rm  "$output"_"$line".mgf; fi; done < "$output"_masses.txt

while read line; do if [ -f  "$output"_"$line"_reconst.txt ]; then rm  "$output"_"$line"_reconst.txt; fi; done < "$output"_masses.txt

rm $allreconstfile
while read line; do if [ -f  "$output"_"$line"_reconst_pvals.txt ]; then cat "$output"_"$line"_reconst_pvals.txt >> $allreconstfile ;rm  "$output"_"$line"_reconst_pvals.txt; fi; done < "$output"_masses.txt

rm "$output"_masses.txt