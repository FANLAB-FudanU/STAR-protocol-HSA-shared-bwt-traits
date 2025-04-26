#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Usage: 
script=scripts/bin/summarize_fm_for_loci_file.sh
bash $script "$per_locus_fm_dir" "$loci_file" "$out_dir" 
'
if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################
per_locus_fm_dir=${1}
loci_file=${2}
out_dir=${3}

echo "per_locus_fm_dir is $per_locus_fm_dir"
echo "loci_file is $loci_file"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

## get chrom_idx start_idx stop_idx
return_ld_prefix=False
read -r chrom_idx start_idx stop_idx <<< "$(get_locus_info_idxs "$loci_file" $return_ld_prefix)"

loginfo "Treating column $chrom_idx as chrom of each locus"
loginfo "Treating column $start_idx as start position of each locus"
loginfo "Treating column $stop_idx as stop position of each locus"

## set output files
mkdir -p "${out_dir}"

loci_with_fm_snps=${out_dir}/loci_with_fm_snps_stats.tsv
echo -e "chrom\tstart\tstop\tnum_CSs\tnum_FM_SNPs\tCS_sizes" > "$loci_with_fm_snps"

cs_sizes_file=${out_dir}/cs_sizes.tsv
echo -e "chrom\tstart\tstop\tcs_id\tcs_size" > "$cs_sizes_file"

num_total_loci=$(sed '1d' "$loci_file" | wc -l)
counter=0
p_iter=100

## get loci with fm snps ans fm stats
loginfo "Getting loci with FM snps and FM stats"
sed '1d' "$loci_file" | awk -vOFS='\t' -v chrom_idx="$chrom_idx" -v start_idx="$start_idx" -v stop_idx="$stop_idx" \
    '{print $chrom_idx, $start_idx, $stop_idx}' | while read -r chrom locus_start locus_stop
do
    locus_id="${chrom}_${locus_start}_${locus_stop}"
    locus_out_dir=${per_locus_fm_dir}/${locus_id}; mkdir -p "${locus_out_dir}"

    locus_fm_out_dir=${locus_out_dir}/outputs
    cred_file=${locus_fm_out_dir}/locus.susie.cred

    if [[ -f $cred_file ]]; then
        fm_snps_file=${locus_fm_out_dir}/locus.susie.snp_in_cs.tsv
        num_CSs=$(sed '1d' "$cred_file" | wc -l)
        num_fm_snps=$(sed '1d' "$fm_snps_file" | wc -l)
        cs_sizes=$(sed '1d' "$cred_file" | awk '{printf("%s,", $6)}')
        echo -e "$chrom\t$locus_start\t$locus_stop\t$num_CSs\t$num_fm_snps\t$cs_sizes" >> "$loci_with_fm_snps"

        # report cs size
        cut -f 1,6 "$cred_file" |sed '1d' |awk -vOFS='\t' -v chrom="$chrom" -v locus_start="$locus_start" -v locus_stop="$locus_stop" '{print chrom, locus_start, locus_stop, $1, $2}' >> "$cs_sizes_file"
    fi

    counter=$((counter+1))

    if [[ $(echo $counter $p_iter | awk '{print $1 % $2}') -eq 0 ]]; then
        loginfo "$counter / $num_total_loci done"
    fi

done 

num_loci_with_fm_snps=$(sed '1d' "$loci_with_fm_snps" |wc -l)

## get stats of credible sets
cs_stats_prefix=${out_dir}/cs_stats
script=${SCRIPT_DIR}/get_fm_cs_stats.R
Rscript "$script" --fm_stats_per_locus_file "$loci_with_fm_snps" --out_prefix "$cs_stats_prefix"

## done
loginfo "Done. loci_with_fm_snps is $loci_with_fm_snps. Total $num_loci_with_fm_snps loci have fine-mapped SNPs"
loginfo "cs_sizes_file is $cs_sizes_file"
loginfo "cs_stats_prefix is $cs_stats_prefix"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
