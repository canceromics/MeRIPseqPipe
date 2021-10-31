#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
#$4 argv 4 : peakCalling_tools_count
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3
peakCalling_tools_count=$4

# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9 
done

function SortTransferBed()
{
    bed_file=$1
    bed_anno_file=$2
    outdir=$3
    ## sort bed by pvalue for rank merge && transfer the origin region of peaks into the bedtools merged region of peaks 
    awk '{ print $1":"$2"-"$3,$5}' ${bed_file} | sort -k1,1 |join -a1 - ${bed_anno_file} | sort -k2,2 -n -r | awk '{print $3}' > ${outdir}/tmp.${bed_file}.location
}
function mergebedByRank()
{
    prefix_id=$1
    out_prefix=$2
    mkdir tmp.${out_prefix}
    cat *_${prefix_id}_*.bed | awk '{print $1"\t"$2*1"\t"$3*1"\t"$1":"$2"-"$3}' > tmp.${out_prefix}/bedtools_${prefix_id}_all_peaks
    sortBed -i tmp.${out_prefix}/bedtools_${prefix_id}_all_peaks |mergeBed -i - -c 4,4 -o collapse,count | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}'  > tmp.${out_prefix}/bedtools_${prefix_id}
    awk -F "\t" '{print $4,$5}' tmp.${out_prefix}/bedtools_${prefix_id} | awk -F '[," "]+' '{for (i=2 ;i<=NF;i++) printf $i" "$1"\n" }' | sort -k1 | uniq > tmp.${out_prefix}/bed_anno_file
    for bedfile in *_${prefix_id}_*.bed
    do
        SortTransferBed $bedfile tmp.${out_prefix}/bed_anno_file tmp.${out_prefix}
    done
    paste -d "\t" tmp.${out_prefix}/tmp*location > ${out_prefix}.bedlist
    peak_number=$(wc -l tmp.${out_prefix}/bed_anno_file | cut -d " " -f 1)
    Rscript merge_peaks_by_rank.R ${out_prefix}.bedlist ${peak_number} ${out_prefix}.bed
    rm -rf tmp.${out_prefix} ${out_prefix}.bedlist
}

# Before merging peaks, normalize all peaks of different tools
for bedfile in *.bed
do
read -u 9
{
    mv $bedfile tmp.$bedfile
    python normalize_peaks.py tmp.$bedfile $bedfile
    rm tmp.$bedfile
    echo >&9
}&
done
wait
# if the number of peakcalling tools > 2
if [ $flag_peakCallingbygroup -gt 0 ]; then
    echo "Start to merge different tools' result of every group"
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebedByRank ${group_id} rank_merged_group_${group_id}
        else
            awk '{OFS="\t";$5=10^-$5;print }' *${group_id}*.bed |sortBed -i - > rank_merged_group_${group_id}.bed
        fi
        echo >&9
    }&
    done
    wait
    mergebedByRank merged_group rank_merged_allpeaks
else
    echo "Start to merge different tools' result of every sample"
    sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_group_id in ${sampleinfo_list}
    do
    read -u 9
    {  
        sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        ## Adding the information of group
        for samplefile in *_${sample_id}_normalized.bed
        do
            mv $samplefile ${samplefile/_normalized.bed/}_${group_id}_normalized.bed
        done
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebedByRank ${sample_id} rank_merged_sample_${group_id}_${sample_id}
        else
            awk '{OFS="\t";$5=10^-$5;print }' *${sample_id}*.bed |sortBed -i - > rank_merged_sample_${group_id}_${sample_id}.bed
        fi
        echo >&9
    }&
    done
    wait
    echo "Start to merge different samples' result of every group"
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        mergebedByRank merged_sample_${group_id} rank_merged_group_${group_id}
        echo >&9
    }&
    done
    wait
    #mergebedByRank merged_sample rank_merged_allpeaks
    cat *_merged_sample_*.bed | sortBed -i - |mergeBed -i - -c 4,5 -o count,mean | awk '$4>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5}' > rank_merged_allpeaks.bed
fi
echo "Rank merged peaks done"


