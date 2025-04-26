#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: "
    echo 'bash $script "$alkws_gwas_file" "$min_info" "$out_file" '
    exit
fi

########################## Parse arguments ##########################
alkws_gwas_file=${1}
min_info=${2}
min_maf=${3}
out_file=${4}

echo "alkws_gwas_file is $alkws_gwas_file"
echo "min_info is $min_info"
echo "min_maf is $min_maf"
echo "out_file is $out_file"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

script=${SCRIPT_DIR}/make_lava_format_sumstats.awk
bgzip -dc "$alkws_gwas_file" | sed '1d' | awk -f "$script" -v min_info="$min_info" -v min_maf="$min_maf" | \
    bgzip > "$out_file"
tabix -S1 -s2 -b3 -e3 "$out_file"

loginfo "Done. out_file is $out_file"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
