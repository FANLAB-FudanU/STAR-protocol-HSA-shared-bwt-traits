#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: "
    echo 'bash $script "$work_dir" "$rsid_file" "$dbsnp_bigbed_file" "$chroms_list" "$out_dir" "$num_rows" "$threads" '
    exit
fi

########################## Parse arguments ##########################
work_dir=${1}
rsid_file=${2}
dbsnp_bigbed_file=${3}
chroms_list=${4}
out_dir=${5}
num_rows=${6}
threads=${7}

echo "work_dir is $work_dir"
echo "rsid_file is $rsid_file"
echo "dbsnp_bigbed_file is $dbsnp_bigbed_file"
echo "chroms_list is $chroms_list"
echo "out_dir is $out_dir"
echo "num_rows is $num_rows"
echo "threads is $threads"

check-files.sh "$work_dir" "$rsid_file" "$dbsnp_bigbed_file" "$chroms_list"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"


mkdir -p "${out_dir}"

### def vars
working_dir=${out_dir}/working; mkdir -p "${working_dir}"

splitted_rsid_dir=${working_dir}/01.splitted_rsid; mkdir -p "${splitted_rsid_dir}"
fetched_dbsnp_dir=${working_dir}/02.fetched_dbsnp; mkdir -p "${fetched_dbsnp_dir}"

fetch_cmds_script=${working_dir}/fetch_cmds.sh
splitted_rsid_prefix=${splitted_rsid_dir}/rsid.

out_dbsnp_file=${out_dir}/dbsnp.bed.gz
out_pos_rsid_ref_allele_file=${out_dir}/snp_pos.rsid.ref_allele.tsv
out_snp_pos_with_chr=${out_dir}/snp_pos.with_chr.tsv
out_snp_pos_without_chr=${out_dir}/snp_pos.without_chr.tsv

### split snp pos into batches
loginfo "Splitting rsid"
split "$rsid_file" "$splitted_rsid_prefix" -l "$num_rows" -d --additional-suffix=.tsv

### gen cmds to fetch dbsnp
loginfo "Generating cmds to fetch dbsnp"

[[ -f $fetch_cmds_script ]] && rm -f "$fetch_cmds_script"
for batch_rsid in "${splitted_rsid_prefix}"*
do
    filename=$(basename "$batch_rsid")
    batch_id=$(echo "$filename" | cut -d. -f2)
    # echo "$batch_id"
    batch_out_file=${fetched_dbsnp_dir}/rsid.${batch_id}.bed
    echo "cd $work_dir" >> "$fetch_cmds_script"
    echo "bigBedNamedItems -nameFile $dbsnp_bigbed_file $batch_rsid $batch_out_file" >> "$fetch_cmds_script"
    echo "grep -Fwf $chroms_list $batch_out_file |cat | bgzip > $batch_out_file.gz" >> "$fetch_cmds_script"
done

### submit jobs
loginfo "Running script ($fetch_cmds_script) to fetch dbsnp"
paralleltask -t local -l 3 -p 1 -m "$threads" --job_prefix fetch_dbsnp --disable_convert_path "$fetch_cmds_script"

### merge fetched rsids
for dbsnp_file in "${fetched_dbsnp_dir}"/rsid*bed.gz
do
    gzip -dc "$dbsnp_file"
done | sort -k1,1 -k2,2n | bgzip > "$out_dbsnp_file"

tabix -s1 -b2 -e3 -0 "$out_dbsnp_file"

##### Get snp pos (1-based), ref allele and rsid from dbsnp155 hg38 bed
loginfo "Getting snp pos, rsid and ref allele"
(
    echo -e "chrom\tpos\tref\trsid"; 
    bgzip -dc "$out_dbsnp_file" | awk -vOFS='\t' '{print $1,$3,$5,$4}'
) > "$out_pos_rsid_ref_allele_file"

sed '1d' "$out_pos_rsid_ref_allele_file" |awk -vOFS='\t' '{print $1, $2, $2}' > "$out_snp_pos_with_chr"
sed 's/chr//' "$out_snp_pos_with_chr" > "$out_snp_pos_without_chr"

### compress cmds and working dir
loginfo "Compressing scripts"
tar -zcf "${out_dir}"/fetch_cmds.tar.gz "${fetch_cmds_script}"*
tar -zcf "${working_dir}".tar.gz "${working_dir}"
rm -rf "$working_dir"

loginfo "Done. out_dbsnp_file is $out_dbsnp_file"
loginfo "out_snp_pos_with_chr is $out_snp_pos_with_chr"
loginfo "out_snp_pos_without_chr is $out_snp_pos_without_chr"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
