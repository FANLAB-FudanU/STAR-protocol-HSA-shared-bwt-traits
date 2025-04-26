#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: "
    echo 'bash $script "$per_locus_fm_dir" "$loci_with_fm_snps" "$out_dir" '
    exit
fi

########################## Parse arguments ##########################
per_locus_fm_dir=${1}
loci_with_fm_snps=${2}
out_dir=${3}

echo "per_locus_fm_dir is $per_locus_fm_dir"
echo "loci_with_fm_snps is $loci_with_fm_snps"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

## get chrom_idx start_idx stop_idx
return_ld_prefix=False
read -r chrom_idx start_idx stop_idx <<< "$(get_locus_info_idxs "$loci_with_fm_snps" $return_ld_prefix)"

loginfo "Treating column $chrom_idx as chrom of each locus"
loginfo "Treating column $start_idx as start position of each locus"
loginfo "Treating column $stop_idx as stop position of each locus"

## def vars
mkdir -p "${out_dir}"
out_merged_snps_file=${out_dir}/merged_FM_SNPs.tsv

## get first fm snps file
read -r chrom locus_start locus_stop <<< "$(sed '1d' "$loci_with_fm_snps" | head -1 | awk -vOFS='\t' -v chrom_idx="$chrom_idx" -v start_idx="$start_idx" -v stop_idx="$stop_idx" '{print $chrom_idx, $start_idx, $stop_idx}')"
locus_id="${chrom}_${locus_start}_${locus_stop}"
locus_out_dir=${per_locus_fm_dir}/${locus_id}
locus_fm_out_dir=${locus_out_dir}/outputs
first_fm_snps=${locus_fm_out_dir}/locus.susie.snp_in_cs.tsv

# get header
head -1 "$first_fm_snps" |awk '{print "chrom\tlocus_start\tlocus_stop\t"$0}' > "$out_merged_snps_file"

## merge all loci
counter=0
p_iter=100
num_total_loci=$(sed '1d' "$loci_with_fm_snps" |wc -l)

loginfo "Merging all loci"
sed '1d' "$loci_with_fm_snps" | awk -vOFS='\t' -v chrom_idx="$chrom_idx" -v start_idx="$start_idx" -v stop_idx="$stop_idx" '{print $chrom_idx, $start_idx, $stop_idx}' | while read -r chrom locus_start locus_stop
do
    locus_id="${chrom}_${locus_start}_${locus_stop}"
    locus_out_dir=${per_locus_fm_dir}/${locus_id}
    locus_fm_out_dir=${locus_out_dir}/outputs

    fm_snps=${locus_fm_out_dir}/locus.susie.snp_in_cs.tsv

    # merge
    sed '1d' "$fm_snps" | awk -v chrom="$chrom" -v locus_start="$locus_start" -v locus_stop="$locus_stop" '{print chrom"\t"locus_start"\t"locus_stop"\t"$0}' >> "$out_merged_snps_file"

    # counter ++
    counter=$((counter+1))
    print_status=$((counter % p_iter))
    if [[ print_status -eq 0 ]]; then
        echo "$counter / $num_total_loci loci are done"
    fi
done

bgzip -f "$out_merged_snps_file"
out_merged_snps_file=${out_merged_snps_file}.gz
tabix -S1 -s5 -b6 -e6 "${out_merged_snps_file}"

loginfo "Getting snp pos"
snp_pos_dir=${out_dir}/snp_pos; mkdir -p "${snp_pos_dir}"
snp_pos_with_chr=${snp_pos_dir}/merged_snps.snp_pos.with_chr.tsv
snp_pos_without_chr=${snp_pos_dir}/merged_snps.snp_pos.without_chr.tsv

bgzip -dc "$out_merged_snps_file" | sed '1d' |awk -vOFS='\t' '{print "chr"$5,$6,$6}' > "$snp_pos_with_chr"
bgzip -dc "$out_merged_snps_file" | sed '1d' |awk -vOFS='\t' '{print $5,$6,$6}' > "$snp_pos_without_chr"

loginfo "Getting avinput"
snpid_avinput_dir=${out_dir}/snp_ids_avinput; mkdir -p "${snpid_avinput_dir}"
avinput=${snpid_avinput_dir}/merged_FM_SNPs.avinput
bgzip -dc "$out_merged_snps_file" | sed '1d' |awk -vOFS='\t' '{print "chr"$5,$6,$6,$7,$8}' > "$avinput"

loginfo "Getting snpid and pos"
snpid_pos_file=${snpid_avinput_dir}/merged_FM_SNPs.ids_pos.tsv
bgzip -dc "$out_merged_snps_file" | cut -f4-6 > "$snpid_pos_file"

loginfo "Done. out_merged_snps_file is ${out_merged_snps_file}. snp_pos_without_chr is $snp_pos_without_chr. "
loginfo "avinput is $avinput"
loginfo "snpid_pos_file is $snpid_pos_file"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
