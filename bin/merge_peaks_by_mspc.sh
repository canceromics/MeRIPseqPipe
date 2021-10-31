#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
#$4 argv 4 : peakCalling_tools_count
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3
peakCalling_tools_count=$4
out_dir=$5
# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9 
done
# Generate the diroectory of results:
mkdir ${out_dir}

# Define the function of MSPC runing for different situations
function mergebedForBio()
{
    prefix_id=$1
    out_prefix=$2
    bedfile_array=$(ls *_${prefix_id}_*.bed | awk '{ORS=" "}{print "-i",$0}')
    mspc -i  $bedfile_array -r bio -s 1E-4 -w 1E-2 -o Bio_$prefix_id
    ln Bio_$prefix_id/ConsensusPeaks.bed ${out_prefix}.bed
    awk 'NR>1{OFS="\t";$5=10^-$5;print $1,$2,$3,$1":"$2"-"$3,$5}' Bio_$prefix_id/ConsensusPeaks.bed |sortBed -i - > ${out_dir}/${out_prefix}.bed
}
function mergebedForTec()
{
    prefix_id=$1
    out_prefix=$2
    peakCalling_tools_count=$3
    bedfile_array=$(ls *_${prefix_id}_*.bed | awk '{ORS=" "}{print "-i",$0}')
    mspc -i $bedfile_array -r tec -s 1E-2 -w 1E-1 -o Tec_$prefix_id
    ln Tec_$prefix_id/ConsensusPeaks.bed ${out_prefix}.bed
    awk 'NR>1{OFS="\t";$5=10^-$5;print $1,$2,$3,$1":"$2"-"$3,$5}' Tec_$prefix_id/ConsensusPeaks.bed |sortBed -i - > ${out_dir}/${out_prefix}.bed
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
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebedForTec ${group_id} mspc_merged_group_${group_id}
        else
            ln *${group_id}*.bed mspc_merged_group_${group_id}.bed
            awk '{OFS="\t";$5=10^-$5;print }' *${group_id}*.bed |sortBed -i - > ${out_dir}/mspc_merged_group_${group_id}.bed
        fi
        echo >&9
    }&
    done
    wait
    mergebedForBio merged_group mspc_merged_allpeaks
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
            mv $samplefile ${samplefile/_normalized.bed/}_${group_id}_normalized.bed
        done
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebedForTec ${sample_id} mspc_merged_sample_${group_id}_${sample_id} $peakCalling_tools_count
        else
            ln *${sample_id}*.bed mspc_merged_sample_${group_id}_${sample_id}.bed
            awk '{OFS="\t";$5=10^-$5;print }' *_${sample_id}_*normalized.bed |sortBed -i - > ${out_dir}/mspc_merged_sample_${group_id}_${sample_id}.bed
        fi
        echo >&9
    }&
    done
    wait
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        mergebedForBio merged_sample_${group_id} mspc_merged_group_${group_id}
        echo >&9
    }&
    done
    wait
    #mergebedForBio merged_sample mspc_merged_allpeaks
    cat ${out_dir}/*_merged_sample_*.bed | sortBed -i - |mergeBed -i - -c 4,5 -o count,mean | awk '$4>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5}' > ${out_dir}/mspc_merged_allpeaks.bed
fi
judge_chr=$(cat *.bed |cut -f 1 |sort |uniq| awk '$0~"chr"{print "includeChr"}' |uniq)
if [ "$judge_chr" != "includeChr" ]; then sed -i 's/chr//g' ${out_dir}/*.bed ;fi
mv ${out_dir}/*.bed ./
exec 9<>-
echo "MSPC merged peaks done"
