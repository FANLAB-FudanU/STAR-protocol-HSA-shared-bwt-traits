#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Usage:
script=scripts/bin/extract_EUR_genotype.sh
bash $script --KGP_Phase3_vcf_dir="$KGP_Phase3_vcf_dir" --sample_info="$sample_info" --threads="$threads" --out_dir="$out_dir" 
'
if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################

while [ $# -gt 0 ]; do
    case "$1" in
	--KGP_Phase3_vcf_dir=*)
		KGP_Phase3_vcf_dir="${1#*=}"
		;;
	--sample_info=*)
		sample_info="${1#*=}"
		;;
	--threads=*)
		threads="${1#*=}"
		;;
	--out_dir=*)
		out_dir="${1#*=}"
		;;
        *)
        echo "* Error: Invalid argument: $1*"
        exit 1
    esac
    shift
done

echo "KGP_Phase3_vcf_dir is $KGP_Phase3_vcf_dir"
echo "sample_info is $sample_info"
echo "threads is $threads"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

function tmp()
{
    KGP_Phase3_vcf_dir=
    sample_info=
    threads=
    out_dir=

}

log_dir=${out_dir}/logs; mkdir -p "${log_dir}"
per_chrom_dir=${out_dir}/per_chrom; mkdir -p "${per_chrom_dir}"

# get EUR sample ids
EUR_sample_ids=${out_dir}/EUR_sample_ids.txt
sed '1d' "$sample_info" |awk '$3=="EUR"{print $1}' |sort -u > "$EUR_sample_ids"

# extract genotype of EUR samples of all chromosomes
sample_ids=${out_dir}/EUR_sample_ids.txt

loginfo "Extracting genotype of EUR samples of all chromosomes"
find "$KGP_Phase3_vcf_dir"/*gz | while read -r vcf_file
do
    chrom=$(echo "$vcf_file" |awk -F'/' '{print $NF}' |awk -F'.' '{print $2}')
    out_prefix=${per_chrom_dir}/$chrom
    log=${log_dir}/extract_$chrom.log
    plink2 --vcf "$vcf_file" --keep "$sample_ids" --snps-only just-acgt \
        --set-missing-var-ids @_#_\$r_\$a --max-alleles 2 --make-pgen vzs --threads "$threads" \
        --out "$out_prefix" > "$log"
    echo "done for $chrom"
done

# merge genotype of EUR samples of all chromosomes
chrom_list=${out_dir}/chrom_list.txt
(
    seq 1 22 |awk '{print "chr"$1}'
    echo "chrX"
) > "$chrom_list"

file_list=${per_chrom_dir}/file_list.txt
[[ -f $file_list ]] && rm -f "$file_list"

while read -r chrom
do
    per_chrom_pfile=${per_chrom_dir}/${chrom}
    echo "$per_chrom_pfile" >> "$file_list"
done < "$chrom_list"

log=${log_dir}/merge_all_geno.log
loginfo "Merging all geno (log: $log)"
out_pfile=${out_dir}/EUR
if ! plink2 --pmerge-list "$file_list" pfile-vzs --threads "$threads" --out "${out_pfile}" \
    --pmerge-output-vzs > "$log" 2>&1; then
    logerror "Error when merging all geno (log: $log)"
fi

# convert to bfile
loginfo "convert to bfile"
plink2 --pfile "${out_dir}"/EUR vzs --make-bed --out "${out_dir}"/EUR --threads "$threads"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
