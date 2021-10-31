#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
#$4 argv 4 : peakCalling_tools_count
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3
peakCalling_tools_count=$4
peakCalling_tools_main=$5

# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9 
done

function mergebedByBedtools()
{
    prefix_id=$1
    out_prefix=$2
    peakCalling_tools_main=$3
    cat ${peakCalling_tools_main}*${prefix_id}*normalized.bed | sortBed -i - | mergeBed -i - -c 4,5 -o count,mean > tmp.${prefix_id}_allPeaks.bed
    ls *${prefix_id}*normalized.bed | grep -v ${peakCalling_tools_main} | xargs -i cat {} | sortBed -i - | mergeBed -i - -c 4,5 -o count,mean > tmp.${prefix_id}_others_allPeaks.bed
    intersectBed -a tmp.${prefix_id}_allPeaks.bed -b tmp.${prefix_id}_others_allPeaks.bed -u | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$5}' > ${out_prefix}.bed
}

if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebedByBedtools ${group_id} bedtools_merged_group_${group_id} ${peakCalling_tools_main}
        else
            awk '{OFS="\t";$5=10^-$5;print }' *${group_id}*.bed | sortBed -i - > bedtools_merged_group_${group_id}.bed
        fi
        echo >&9
    }&
    done
    wait
    if [ $peakCalling_tools_count -gt 1 ]; then
        mergebedByBedtools "" bedtools_merged_allpeaks ${peakCalling_tools_main}
    else
        cat ${peakCalling_tools_main}_*_normalized.bed | sortBed -i - | mergeBed -i - -c 4,5 -o count,mean | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$5}' > bedtools_merged_allpeaks.bed
    fi
else
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
            mv $samplefile ${samplefile/_normalized.bed/_${group_id}_normalized.bed}
        done
        echo >&9
    }&
    done
    wait
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebedByBedtools ${group_id} bedtools_merged_group_${group_id} ${peakCalling_tools_main}
        else
            awk '{OFS="\t";$5=10^-$5;print }' *${group_id}*.bed | sortBed -i - > bedtools_merged_group_${group_id}.bed
        fi
        echo >&9
    }&
    done
    wait
    if [ $peakCalling_tools_count -gt 1 ]; then
        mergebedByBedtools "" bedtools_merged_allpeaks ${peakCalling_tools_main}
    else
        cat ${peakCalling_tools_main}_*_normalized.bed | sortBed -i - | mergeBed -i - -c 4,5 -o count,mean | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$5}' > bedtools_merged_allpeaks.bed
    fi
fi
echo "${peakCalling_tools_main} merged peaks done"
