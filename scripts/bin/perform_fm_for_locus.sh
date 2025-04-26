#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Usage:
script=scripts/bin/perform_fm_for_locus.sh
bash $script "$chrom" "$locus_start" "$locus_stop" "$gwas_sumstats" "$ld_prefixes_file" \
    "$min_freq" "$min_info" "$max_pval" "$n_samples" "$n_causal_snps" "$out_dir" 

'

if [ $# -eq 0 ]; then
    echo "Usage: "
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################
chrom=${1}
locus_start=${2}
locus_stop=${3}
gwas_sumstats=${4}
ld_prefixes_file=${5}
min_freq=${6}
min_info=${7}
max_pval=${8}
n_samples=${9}
n_causal_snps=${10}
out_dir=${11}

echo "chrom is $chrom"
echo "locus_start is $locus_start"
echo "locus_stop is $locus_stop"
echo "gwas_sumstats is $gwas_sumstats"
echo "ld_prefixes_file is $ld_prefixes_file"
echo "min_freq is $min_freq"
echo "min_info is $min_info"
echo "max_pval is $max_pval"
echo "n_samples is $n_samples"
echo "n_causal_snps is $n_causal_snps"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

mkdir -p "${out_dir}"
fm_inputs_dir=${out_dir}/inputs; mkdir -p "${fm_inputs_dir}"
fm_outputs_dir=${out_dir}/outputs; mkdir -p "${fm_outputs_dir}"

loginfo "Preparing fine-mapping inputs"

script=${SCRIPT_DIR}/prepare_fm_inputs_for_locus.py
python -u "$script" --chrom "$chrom" \
    --locus-start "$locus_start" --locus-stop "$locus_stop" \
    --gwas-smmr "$gwas_sumstats" --ld-prefixes-file "$ld_prefixes_file" \
    --min-freq "$min_freq" --min-info "$min_info" --max-pval "$max_pval" \
    --out-dir "$fm_inputs_dir"

## run susie
zfile=${fm_inputs_dir}/03.z_ld/sumstats.z
ldfile=${fm_inputs_dir}/03.z_ld/sumstats.ld

fm_outputs_prefix=${fm_outputs_dir}/locus

loginfo "Running SuSiE. fm_outputs_prefix is $fm_outputs_prefix"
script=${SCRIPT_DIR}/run_susieR.R
Rscript "$script" --z "$zfile" --ld "$ldfile" \
    --out "$fm_outputs_prefix" --n-samples "$n_samples" --L "$n_causal_snps" \
    --save-susie-obj 

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"






