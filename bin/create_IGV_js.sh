#! /bin/bash
fasta=$1
gtf=$2
merged_peak_file=$3
designfile=$4
echo "Start to generate IGV.js"

## setting tmp files' name
bedgraph_tracks_file=tmp.bedgraph.tracks
peaks_tracks_file=tmp.peaks.tracks

## combined tracks of bedgraph
sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
for sample_group_id in ${sampleinfo_list}
    do
    {  
        sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        bedgraph_input_file=$(ls ${sample_id}.input_*.igv.bedgraph)
        bedgraph_ip_file=$(ls ${sample_id}.ip_*.igv.bedgraph)
        cat >> ${bedgraph_tracks_file} << EOF
        {
            url: '${bedgraph_input_file}',
            name: '${sample_id}.input',
            color: 'rgb(200,0,0)',
            type: "wig",
            sourceType: "file",
            autoscaleGroup: 'group_${group_id}.${sample_id}'
        },
        {
            url: '${bedgraph_ip_file}',
            name: '${sample_id}.ip',
            type: "wig",
            sourceType: "file",
            color: 'rgb(200,0,0)',
            autoscaleGroup: 'group_${group_id}.${sample_id}'
        },
EOF
    }
done

## combined tracks of merged group peaks
groups_peak_file=$(ls *_merged_group_*igv.bed)
for peak_file in ${groups_peak_file}
    do
    {  
        cat >> ${peaks_tracks_file} << EOF
        {
            type: "annotation",
            format: "bed",
            url: '${peak_file}',
            name: "${peak_file}"
        },
EOF
    }
done

## combined tracks and allpeaks track
cat ${bedgraph_tracks_file} ${peaks_tracks_file} > tmp.tracks
cat >> tmp.tracks << EOF
        {
            type: "annotation",
            format: "bed",
            url: '${merged_peak_file}',
            name: "${merged_peak_file}"
        }
EOF
tracks_js=$(cat tmp.tracks)

## combined all info 
cat>igv.js<<EOF
var igvDiv = document.getElementById("igvDiv");
var options =
{
    reference: {
        id: "$fasta",
        fastaURL: "$fasta",
        indexURL: "$fasta.fai",
        wholeGenomeView: false
    },
    locus: 'chr22',
    tracks: [
        {
            type: "annotation",
            format: "gtf",
            url: "$gtf",
            displayMode: "SQUISHED",
            name: "$gtf",
            visibilityWindow: 10000000
        },
$tracks_js
]
}; 
var browser = igv.createBrowser(igvDiv, options);
EOF
rm -rf tmp*
echo "Generate IGV.js was success"
