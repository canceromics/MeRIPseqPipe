#!/bin/bash
#bash m6Aprediction.sh <matk_jar> <fasta> <gtf> <THREAD_NUM>
## $1 argv 1 : matk_jar
## $2 argv 2 : designfile
## $3 argv 3 : fasta file
## $4 argv 4 : gtf file
matk_jar=$1
designfile=$2
fasta_file=$3
gtf_file=$4

### check if the file matk.jar exists
if [ ! -f "$matk_jar" ]; then
    echo "Cannot find matk.jar. Please check the param of matk_jar" 1>&2
    exit 1
fi

faToTwoBit ${fasta_file} ${fasta_file/.fa/.2bit}
awk -F "\t" '$3=="gene"{print }' $gtf_file > tmp.$gtf_file
group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
for group_id in $group_list
do 
{
    bedfile=$(ls *merged_group_${group_id}.bed)
    ip_bam_file_array=$(echo *.ip_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    input_bam_file_array=$(echo *.input_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    java -jar $matk_jar -singleNucleotide \
                    -mode "MeRIP" \
                    -ip "$ip_bam_file_array" \
                    -input "$input_bam_file_array" \
                    -bed ${bedfile} \
                    -2bit ${fasta_file/.fa/.2bit} \
                    -gtf tmp.$gtf_file \
                    -out m6A_sites_${group_id}.bed
    awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t*\t"$6"\t"$5"\t'${group_id}'"}' m6A_sites_${group_id}.bed > tmp.m6A_sites_${group_id}.bed
}
done
wait
cat tmp.m6A_sites*.bed | sortBed | mergeBed -s -c 4,6,7,8 -o first,first,collapse,collapse > tmp.m6A_sites_merged.bed
awk -v gap=25 '{print $1"\t"$2-gap"\t"$3+gap"\t*\t*\t"$5}' tmp.m6A_sites_merged.bed | bedtools getfasta -s -fi ${fasta_file} -bed - | awk '$0!~">"{print $0}' > tmp.m6A_sites_merged.fa
paste tmp.m6A_sites_merged.bed tmp.m6A_sites_merged.fa > m6A_sites_merged.bed
rm tmp.*
echo "Prediction sites of m6A done"
