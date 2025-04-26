#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Main: Extract MAFs for SNPs and reformat MAFs. 
Usage: 
Args:
    raw_maf_dir: Directory contains raw maf files. Default is /public/home/fan_lab/wangjie/PublicDB/UCSC/01.DownloadedData/multiz30way_hg38/
    snp_pos: 0-based position of SNPs (bed file). No header. Three columns. Chrom column contains chr string. Only first three columns are used. Additional columns are allowed but not used
    human_assembly: hg38. 

bash $script "$work_dir" "$raw_maf_dir" "$snp_pos" "$human_assembly" "$species_assembly_file" "$threads" "$out_prefix" 
'

if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################
work_dir=${1}
raw_maf_dir=${2}
snp_pos=${3}
human_assembly=${4}
species_assembly_file=${5}
threads=${6}
out_prefix=${7}

echo "work_dir is $work_dir"
echo "raw_maf_dir is $raw_maf_dir"
echo "snp_pos is $snp_pos"
echo "human_assembly is $human_assembly"
echo "species_assembly_file is $species_assembly_file"
echo "threads is $threads"
echo "out_prefix is $out_prefix"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

# def out files
out_working_dir=${out_prefix}.working; mkdir -p "${out_working_dir}"
splitted_snp_pos_dir=${out_working_dir}/01.splitted_snp_pos; mkdir -p "${splitted_snp_pos_dir}"
fetched_maf_dir=${out_working_dir}/02.fetched_maf; mkdir -p "${fetched_maf_dir}"

# split snp pos by chrom
loginfo "Spliting snp pos by chrom"
chroms_list=${splitted_snp_pos_dir}/chroms_list.txt
cut -f1 "$snp_pos" | sort |uniq > "$chroms_list"

while read -r chrom
do
    echo "chrom is $chrom"
    awk -v chrom="$chrom" '$1==chrom' "$snp_pos" > "${splitted_snp_pos_dir}"/"${chrom}".pos.bed
done < "$chroms_list"

loginfo "Done. splitted_snp_pos_dir is $splitted_snp_pos_dir"

# fetch maf for each chrom
loginfo "Fetching maf for each chrom"

script=$SCRIPT_DIR/fetch_maf_by_snp_pos.py

cmds_script=${out_working_dir}/fetch_maf_for_each_chrom_cmds.sh
[[ -f $cmds_script ]] && rm -rf "$cmds_script"*

cat "$chroms_list" | while read -r chrom
do
    raw_maf_file=${raw_maf_dir}/${chrom}.maf.gz
    snp_pos_file="${splitted_snp_pos_dir}"/"${chrom}".pos.bed
    chrom_maf_prefix=${fetched_maf_dir}/${chrom}

    (
        echo "cd $work_dir"; 
        echo "python -u $script --raw_maf_file $raw_maf_file --snp_pos_file $snp_pos_file \
            --human_assembly $human_assembly --out_prefix $chrom_maf_prefix"
    )
done > "$cmds_script"

loginfo "Executing $cmds_script to fetch maf"

paralleltask -t local -l 2 -p 1 --rerun 0 -m "$threads" --disable_convert_path --job_prefix fetch_maf "$cmds_script"
rm -f pid*log.info

loginfo "MAF are stored in $fetched_maf_dir"

fetched_maf_file=${out_prefix}.maf.tsv.gz
out_status_file=${out_prefix}.snp_status.tsv.gz

# merge mafs of all chrom into one file
loginfo "Merging maf of each chrom into one file"
(
    echo '##maf version=1';
    cat "$chroms_list" | while read -r chrom
    do
        chrom_maf_file=${fetched_maf_dir}/${chrom}.maf.tsv.gz
        gzip -dc "$chrom_maf_file" | sed '1d'
    done 
) | gzip > "${fetched_maf_file}"

# reformat raw mafs
loginfo "Reformatting raw mafs"
reformatted_maf_file=${out_prefix}.refmtted.maf.tsv.gz

script=${SCRIPT_DIR}/reformat_multizalign_file.py
python -u "$script" --input_file "$fetched_maf_file" --species_assembly_file "$species_assembly_file" \
    --human_assembly "$human_assembly" --out_file "$reformatted_maf_file"

loginfo "Merging status file into one file"
# merge status file
(
    echo -e 'chrom\tstart\tend\tfound_maf';
    cat "$chroms_list" | while read -r chrom
    do
        chrom_status_file=${fetched_maf_dir}/${chrom}.snp_status.tsv
        sed '1d' "$chrom_status_file"
    done 
) | gzip > "${out_status_file}"

# 打包working目录
tar -zcf "${out_working_dir}".tar.gz "$out_working_dir"
rm -rf "$out_working_dir"

loginfo "Done. fetched_maf_file is $fetched_maf_file"
loginfo "reformatted_maf_file is $reformatted_maf_file"
loginfo "out_status_file is $out_status_file"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
