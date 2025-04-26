#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Usage:
script=scripts/bin/prepare_hsa_inputs.sh
bash $script --hg19_to_hg38_coord_lookup_table="$hg19_to_hg38_coord_lookup_table" --isec_fm_snps_file="$isec_fm_snps_file" --out_file="$out_file" 
'
if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################

while [ $# -gt 0 ]; do
    case "$1" in
	--hg19_to_hg38_coord_lookup_table=*)
		hg19_to_hg38_coord_lookup_table="${1#*=}"
		;;
	--isec_fm_snps_file=*)
		isec_fm_snps_file="${1#*=}"
		;;
	--out_file=*)
		out_file="${1#*=}"
		;;
        *)
        echo "* Error: Invalid argument: $1*"
        exit 1
    esac
    shift
done

echo "hg19_to_hg38_coord_lookup_table is $hg19_to_hg38_coord_lookup_table"
echo "isec_fm_snps_file is $isec_fm_snps_file"
echo "out_file is $out_file"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

function tmp()
{
    hg19_to_hg38_coord_lookup_table=
    isec_fm_snps_file=
    out_file=

}

(
    echo -e "chrom\tstart\tend\tref\talt";
    bgzip -dc $hg19_to_hg38_coord_lookup_table |sed '1d' | awk -F'\t' -vOFS='\t' 'NR==FNR{
        rsid=$1;
        hg38_pos=$5;
        hg38_ref=$6;
        rsid_to_hg38_pos[rsid]=hg38_pos;
        rsid_to_hg38_ref[rsid]=hg38_ref
    } NR>FNR {
        rsid=$4; chrom=$5; hg19_ref=$7; hg19_alt=$8;
        hg38_pos=rsid_to_hg38_pos[rsid];
        hg38_ref=rsid_to_hg38_ref[rsid];
        if (hg38_ref == hg19_ref) {
            hg38_alt=hg19_alt
        } else {
            hg38_alt=hg19_ref
        }
        print "chr"chrom, hg38_pos, hg38_pos, hg38_ref, hg38_alt
    }' - <(bgzip -dc $isec_fm_snps_file |sed '1d')
) > $out_file

loginfo "Done. out_file is $out_file"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
