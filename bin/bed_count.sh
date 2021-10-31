#!/bin/bash
#bash bed_count.sh <designfile> <THREAD_NUM> <merge_bed_file> <output_bam_stat_file>
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : merge_bed_file
#$4 argv 4 : output_bam_stat_file
designfile=$1
THREAD_NUM=$2
merge_bed_file=$3
output_bam_stat_file=$4

# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

# Create the file about the summary of bam stat
echo "Total_Reads" > $output_bam_stat_file
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ${merge_bed_file} > tmp.${merge_bed_file}

sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
for sample_group_id in ${sampleinfo_list}
do
read -u 9
{
    sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
    group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}') 

    #Define the input/ip file of the sample
    input_bam_file=$(ls ${sample_id}.input*.bam | awk '{ORS=" "}{print $0}')
    ip_bam_file=$(ls ${sample_id}.ip*.bam | awk '{ORS=" "}{print $0}')

    ## Create files and print the name of samples && Get the Total reads of the samples of input and ip
    echo -e ${input_bam_file}"\t" | awk 'BEGIN{ORS=""}{print $0}' > ${sample_id}.bam_stat.txt
    samtools view -c ${input_bam_file} >> ${sample_id}.bam_stat.txt
    echo -e ${ip_bam_file}"\t" | awk 'BEGIN{ORS=""}{print $0}' >> ${sample_id}.bam_stat.txt
    samtools view -c ${ip_bam_file} >> ${sample_id}.bam_stat.txt

    ## Setting colnames of peaks input/ip count
    echo $input_bam_file \
    | awk 'BEGIN{ORS=""}{print "chrom\tchromStart\tchromEND\tPeakName\t"}{for(x=1;x<NF;x++) print $x"\t" }END{print $x"\n"}' \
    > ${merge_bed_file}.${group_id}.${sample_id}.input.count
    echo ${ip_bam_file} \
    | awk 'BEGIN{ORS=""}{print "chrom\tchromStart\tchromEND\tPeakName\t"}{for(x=1;x<NF;x++) print $x"\t" }END{print $x"\n"}' \
    > ${merge_bed_file}.${group_id}.${sample_id}.ip.count

    ## Count input/ip peaks
    
    bedtools multicov -bams ${input_bam_file} -bed tmp.${merge_bed_file} >> ${merge_bed_file}.${group_id}.${sample_id}.input.count
    bedtools multicov -bams ${ip_bam_file} -bed tmp.${merge_bed_file} >> ${merge_bed_file}.${group_id}.${sample_id}.ip.count
    echo >&9

    # awk -v bam="$input_bam" -v pre="$prefix" '
    # {print " bedtools multicov -bams '${input_bam_file}' -bed tmp.'${merge_bed_file}' >> '${merge_bed_file}'.'${group_id}'.'${sample_id}'.input.count; \
    # bedtools multicov -bams '${ip_bam_file}' -bed tmp.'${merge_bed_file}' >> '${merge_bed_file}'.'${group_id}'.'${sample_id}'.ip.count; \
    # sortBed -i ./"pre".tmp/input/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/input/"$1".bin25.txt"}' $chrName_file \
    # | xargs -iCMD -P$THREAD_NUM bash -c CMD
}&
done
wait
cat *.bam_stat.txt >> $output_bam_stat_file
wait
echo "bedtools count done"
