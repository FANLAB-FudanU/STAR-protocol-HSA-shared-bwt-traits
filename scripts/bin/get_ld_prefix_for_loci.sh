#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Usage:
script=scripts/bin/get_ld_prefix_for_loci.sh
bash $script "$ld_dir" "$loci_file" "$out_prefix" 
'
if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################
ld_dir=${1}
loci_file=${2}
out_prefix=${3}

echo "ld_dir is $ld_dir"
echo "loci_file is $loci_file"
echo "out_prefix is $out_prefix"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

loci_file_with_ld_prefix=${out_prefix}.with_ld_status.tsv
loci_file_with_npz_ld_prefix=${out_prefix}.with_npz_ld.tsv

## add LD status to loci file
script=${SCRIPT_DIR}/get_ld_prefix_for_loci.py
python -u "$script" --ld-dir "$ld_dir" --loci-file "$loci_file" \
    --out-file "$loci_file_with_ld_prefix"

## get loci with LD
(
    head -1 "$loci_file_with_ld_prefix";
    sed '1d' "$loci_file_with_ld_prefix" | awk -F'\t' '$10=="npz_file"'
) > "$loci_file_with_npz_ld_prefix"

loginfo "Done. Loci file with LD file status is $loci_file_with_ld_prefix"
loginfo "Loci file with npz LD file is $loci_file_with_npz_ld_prefix"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
