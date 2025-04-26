#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Main: Get genomic positions in hg19 and hg38 versions from given RSID file. 

Args:
    work_dir: working directory
    dbsnp_hg38_bigbed_file: Path to dbSnp bigbed file in hg38 version. Default: /public/home/fan_lab/wangjie/PublicDB/UCSC/01.DownloadedData/hg38/dbsnp/dbSnp155.hg38.bb. 
    dbsnp_hg19_bigbed_file: Path to dbSnp bigbed file in hg19 version. Default: /public/home/fan_lab/wangjie/PublicDB/UCSC/01.DownloadedData/hg19/dbsnp/dbSnp155.hg19.bb. 
    rsid_file: File containing RSIDs per line. 
    chroms_list: File containing chromosomes (with chr). Only SNPs in these chromosomes will be retained. 
    num_rows: Batch size to fetch SNPs from `dbsnp_hg38_bigbed_file` and `dbsnp_hg19_bigbed_file`. Typically 1000 to 10000 would be a good choice.
    threads: Number of threads used when fetching SNPs from `dbsnp_hg38_bigbed_file` and `dbsnp_hg19_bigbed_file`. 
    out_dir: Output directory with three directories: 01.hg19_pos, 02.hg38_pos and 03.coord_lookup_table. 
'

if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################
work_dir=${1}
dbsnp_hg38_bigbed_file=${2}
dbsnp_hg19_bigbed_file=${3}
rsid_file=${4}
chroms_list=${5}
num_rows=${6}
threads=${7}
out_dir=${8}

echo "work_dir is $work_dir"
echo "dbsnp_hg38_bigbed_file is $dbsnp_hg38_bigbed_file"
echo "dbsnp_hg19_bigbed_file is $dbsnp_hg19_bigbed_file"
echo "rsid_file is $rsid_file"
echo "chroms_list is $chroms_list"
echo "num_rows is $num_rows"
echo "threads is $threads"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

mkdir -p "${out_dir}"

hg19_pos_out_dir=${out_dir}/01.hg19_pos; mkdir -p "${hg19_pos_out_dir}"
hg38_pos_out_dir=${out_dir}/02.hg38_pos; mkdir -p "${hg38_pos_out_dir}"

### 1. Get hg19 pos
loginfo "1 | Getting hg19 pos"
script=${SCRIPT_DIR}/get_snp_pos_by_rsid.sh
bash "$script" "$work_dir" "$rsid_file" "$dbsnp_hg19_bigbed_file" "$chroms_list" "$hg19_pos_out_dir" "$num_rows" "$threads"

### 2. Get hg38 pos
loginfo "2 | Getting hg38 pos"
script=${SCRIPT_DIR}/get_snp_pos_by_rsid.sh
bash "$script" "$work_dir" "$rsid_file" "$dbsnp_hg38_bigbed_file" "$chroms_list" "$hg38_pos_out_dir" "$num_rows" "$threads"

### 3. Merge hg19 pos with hg38 pos
# output format: rsid, chrom, hg19_pos, hg19_ref, hg38_pos, hg38_ref
loginfo "3 | Merging hg19 pos with hg38 pos"
hg19_pos_rsid_file=${hg19_pos_out_dir}/snp_pos.rsid.ref_allele.tsv
hg38_pos_rsid_file=${hg38_pos_out_dir}/snp_pos.rsid.ref_allele.tsv

coord_lookup_table_dir=${out_dir}/03.coord_lookup_table; mkdir -p "${coord_lookup_table_dir}"
hg19_to_hg38_coord_lookup_table_unsorted=${coord_lookup_table_dir}/hg19_to_hg38_coord_lookup_table.with_chr.unsorted.tsv
hg19_to_hg38_coord_lookup_table=${coord_lookup_table_dir}/hg19_to_hg38_coord_lookup_table.with_chr.tsv.gz
hg38_to_hg19_coord_lookup_table=${coord_lookup_table_dir}/hg38_to_hg19_coord_lookup_table.with_chr.tsv.gz

(
    set -ue; 
    echo -e "rsid\tchrom\thg19_pos\thg19_ref\thg38_pos\thg38_ref"; 
    awk -F'\t' -vOFS='\t' 'NR==FNR && FNR>1{
        rsid=$4; hg19_pos=$2; hg19_ref=$3; 
        rsid2hg19pos[rsid]=hg19_pos; rsid2hg19ref[rsid]=hg19_ref;
    } NR > FNR && FNR>1 {
        rsid=$4; hg19_ref=rsid2hg19ref[rsid]; hg19_pos=rsid2hg19pos[rsid]; 
        if (hg19_ref ==  "") {
            printf("Warning! SNP %s do not have hg19 pos", rsid); 
            next 
        } else {
            chrom=$1;
            hg38_pos=$2; hg38_ref=$3; 
            printf("%s\t%s\t%s\t%s\t%s\t%s\n", rsid, chrom, hg19_pos, hg19_ref, hg38_pos, hg38_ref)
        }
    }' "$hg19_pos_rsid_file" "$hg38_pos_rsid_file"
)  > "$hg19_to_hg38_coord_lookup_table_unsorted"

(
    head -1 "$hg19_to_hg38_coord_lookup_table_unsorted"; 
    sed '1d' "$hg19_to_hg38_coord_lookup_table_unsorted" |sort -k2,2 -k3,3n 
) | bgzip > "$hg19_to_hg38_coord_lookup_table"
tabix -S1 -s2 -b3 -e3 "$hg19_to_hg38_coord_lookup_table"

(
    bgzip -dc "$hg19_to_hg38_coord_lookup_table" | head -1;
    bgzip -dc "$hg19_to_hg38_coord_lookup_table" | sed '1d' |sort -k2,2 -k5,5n
) | bgzip > "$hg38_to_hg19_coord_lookup_table"
tabix -S1 -s2 -b5 -e5 "$hg38_to_hg19_coord_lookup_table"

loginfo "Done."
loginfo "hg19_pos_out_dir is $hg19_pos_out_dir"
loginfo "hg38_pos_out_dir is $hg38_pos_out_dir"
loginfo "hg19_to_hg38_coord_lookup_table is $hg19_to_hg38_coord_lookup_table"
loginfo "hg38_to_hg19_coord_lookup_table is $hg38_to_hg19_coord_lookup_table"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
