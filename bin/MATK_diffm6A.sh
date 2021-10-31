#!/bin/bash
## MATK_diffpeakCalling.sh <matk_jar> <designfile> <gtf> <compare_str>
## $1 argv 1 : matk_jar
## $2 argv 2 : designfile
## $3 argv 3 : gtf file
## $4 argv 4 : compare_str
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
matk_jar=$1
designfile=$2
gtf_file=$3
compare_str=$4
merged_bed=$5

# setting the function of Running the quantification mode of MATK by two names of groups
function matk_diffm6a_by_two_id()
{
    group_id_1=$1
    group_id_2=$2
    matk_jar=$3
    gtf_file=$4
    control_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    control_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_ip_bam_file_array=$(echo *ip_${group_id_2}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_input_bam_file_array=$(echo *input_${group_id_2}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    java -jar ${matk_jar} -diff \
                    -control_ip "${control_ip_bam_file_array}" \
                    -control_input "${control_input_bam_file_array}" \
                    -treated_ip "${treated_ip_bam_file_array}" \
                    -treated_input "${treated_input_bam_file_array}" \
                    -control_bed $merged_bed \
                    -treated_bed $merged_bed \
                    -gtf ${gtf_file} \
                    -out tmp.${group_id_1}_${group_id_2}.txt
    awk 'BEGIN{print "ID","Gene_symbol","ControlScore","TreatScore","log2FC","pvalue","qvalue"} \
        NR>1{print $1":"$2"-"$3,$4,$5,$6,$7,log($8)/log(2),$9,$10}' tmp.${group_id_1}_${group_id_2}.txt  > MATK_diffm6A_${group_id_1}_${group_id_2}.txt
}

if [ "$compare_str" != "two_group" ]; then
    # Running MATK quantification with compare_str
    group_id_1=$(echo $compare_str | awk 'BEGIN{FS="_vs_"}{print $1}')
    group_id_2=$(echo $compare_str | awk 'BEGIN{FS="_vs_"}{print $2}')
    matk_diffm6a_by_two_id $group_id_1 $group_id_2 $matk_jar $gtf_file
else
    # Running MATK quantification without compare_str beacause of only two groups
    echo "no compare file"
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS="\t"}{print $0}')
    group_id_1=$(echo $group_list | awk 'BEGIN{FS="\t"}{print $1}')
    group_id_2=$(echo $group_list | awk 'BEGIN{FS="\t"}{print $2}')   
    matk_diffm6a_by_two_id $group_id_1 $group_id_2 $matk_jar $gtf_file
fi
wait
echo "diffMATK done"
