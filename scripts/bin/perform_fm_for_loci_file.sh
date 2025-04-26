#!/usr/bin/bash 

## perform fine-mapping for loci from loci file
## requires: scripts/bin/perform_fm_for_locus.sh

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: "
    echo 'script=scripts/bin/perform_fm_for_loci_file.sh'
    echo 'bash $script "$work_dir" "$loci_file" "$gwas_sumstats" "$min_freq" "$min_info" "$max_pval" "$n_samples" "$n_causal_snps" "$threads" "$dry_run" "$out_dir" '
    exit
fi

########################## Parse arguments ##########################
work_dir=${1}
loci_file=${2} # loci file with start and stop coord and ld_prefix added
gwas_sumstats=${3}
min_freq=${4}
min_info=${5}
max_pval=${6}
n_samples=${7}
n_causal_snps=${8}
threads=${9}
dry_run=${10} # True or False
out_dir=${11}

echo "work_dir is $work_dir"
echo "loci_file is $loci_file"
echo "gwas_sumstats is $gwas_sumstats"
echo "min_freq is $min_freq"
echo "min_info is $min_info"
echo "max_pval is $max_pval"
echo "n_samples is $n_samples"
echo "n_causal_snps is $n_causal_snps"
echo "threads is $threads"
echo "dry_run is $dry_run"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

## get chrom_idx start_idx stop_idx ld_prefix_idx
return_ld_prefix=True
read -r chrom_idx start_idx stop_idx ld_prefix_idx <<< "$(get_locus_info_idxs "$loci_file" $return_ld_prefix)"

loginfo "Treating column $chrom_idx as chrom of each locus"
loginfo "Treating column $start_idx as start position of each locus"
loginfo "Treating column $stop_idx as stop position of each locus"
loginfo "Treating column $ld_prefix_idx as ld_prefix of each locus"

## set out dir
mkdir -p "${out_dir}"
per_locus_out_dir=${out_dir}/per_locus; mkdir -p "${per_locus_out_dir}"

## generate cmds
loginfo "Generating cmds to perform finemapping"
cmds_script=${out_dir}/cmds_to_perform_fm.sh
[[ -f $cmds_script ]] && rm -f "$cmds_script"

script=${SCRIPT_DIR}/perform_fm_for_locus.sh

sed '1d' "$loci_file" | awk -vOFS='\t' -v chrom_idx="$chrom_idx" -v start_idx="$start_idx" -v stop_idx="$stop_idx" -v ld_prefix_idx="$ld_prefix_idx" \
    '{print $chrom_idx, $start_idx, $stop_idx, $ld_prefix_idx}' | while read -r chrom locus_start locus_stop ld_prefix
do
    locus_id="${chrom}_${locus_start}_${locus_stop}"
    locus_out_dir=${per_locus_out_dir}/${locus_id}; mkdir -p "${locus_out_dir}"

    ld_prefixes_file=${locus_out_dir}/ld_prefixes.tsv
    echo "$ld_prefix" | sed 's/;/\n/g' > "$ld_prefixes_file"

    echo "cd $work_dir" >> "$cmds_script"
    echo "bash $script $chrom $locus_start $locus_stop $gwas_sumstats \
        $ld_prefixes_file $min_freq $min_info $max_pval \
        $n_samples $n_causal_snps $locus_out_dir " >> "$cmds_script"
done

if [[ $dry_run == "True" ]]; then
    echo "Dry run mode is on"
    echo "Cmds have been written into $cmds_script"
    echo "Bye!"
    exit 1
else
    ## perform finemapping for each locus
    loginfo "Parallel running $cmds_script"
    paralleltask -t local --lines 2 --cpu 1 --rerun 1 --maxjob "$threads" --job_prefix finemapping --disable_convert_path "$cmds_script"
fi

loginfo "Done. out_dir is $out_dir"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
