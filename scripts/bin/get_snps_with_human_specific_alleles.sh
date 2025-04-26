#!/usr/bin/bash 

set -ue
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)
source "${SCRIPT_DIR}/utils.sh"

usage_msg='Main: Get SNPs with human specific alleles. Inputs and outputs are both hg38 version. 

Usage:
    log=tmp/1.log
    script=scripts/human_specific_alleles/get_snps_with_human_specific_alleles.sh
    bash $script "$work_dir" "$human_alleles_file" "$raw_maf_dir" "$human_assembly" "$species_assembly_file" "$threads" "$out_dir" > $log 2>&1 &

Args:
    work_dir: Working directory. 
    human_alleles_file: Five columns: chrom, start, end, ref and alt. Position is 1-based in hg38 version. Chrom column contains chr string. Header should be: chrom, start, end, ref and alt.
    raw_maf_dir: Directory containing MAFs files separated by chromosome in the format of {chrom}.maf.gz. For example, chr1.maf.gz. Default: /public/home/fan_lab/wangjie/PublicDB/UCSC/01.DownloadedData/hg38/multiz30way
    human_assembly: Assembly version of human genome position. Currently only hg38 is supported. 
    species_assembly_file: The species name to assembly version table where assembly version must be contained by MAFs from `raw_maf_dir`. Default: /public/home/fan_lab/wangjie/PublicDB/AncientGenome/040.SNPsAnalysis/030.ExternalSources/010.UCSC_MultizAlign/cons30primates_hg38.non_human_primates.lookup_table.tsv
    threads: Number threads to parallel this process. 
    out_dir: Output directory. Genome build of output files is hg38. 
'

if [ $# -eq 0 ]; then
    echo -e "$usage_msg"
    exit
fi

########################## Parse arguments ##########################
work_dir=${1}
human_alleles_file=${2}
raw_maf_dir=${3}
human_assembly=${4}
species_assembly_file=${5}
threads=${6}
out_dir=${7}

echo "work_dir is $work_dir"
echo "human_alleles_file is $human_alleles_file"
echo "raw_maf_dir is $raw_maf_dir"
echo "human_assembly is $human_assembly"
echo "species_assembly_file is $species_assembly_file"
echo "threads is $threads"
echo "out_dir is $out_dir"

########################## Begin ##########################
starttime=$(date +"%Y-%m-%d %H:%M:%S")
loginfo "Begin to analysis"

mkdir -p "${out_dir}"

### 0. def output vars
snp_maf_dir=${out_dir}/01.maf_of_snps; mkdir -p "${snp_maf_dir}"
snp_pos_bed=${snp_maf_dir}/snps.pos.with_chr.bed

merged_alleles_dir=${out_dir}/02.merged_alleles; mkdir -p "${merged_alleles_dir}"

### 1. get bed pos
less "$human_alleles_file" |sed '1d' |awk -vOFS='\t' '{print $1,$2-1,$3}' > "$snp_pos_bed"

### 2. fetch MAFs for SNPs from raw MAFs and reformat MAF
snp_maf_prefix=${snp_maf_dir}/snps

script=$SCRIPT_DIR/fetch_multizalign_from_maf.sh
bash "$script" "$work_dir" "$raw_maf_dir" "$snp_pos_bed" "$human_assembly" "$species_assembly_file" "$threads" "$snp_maf_prefix" 

### 3. merge reformatted MAF with human alleles
refmtted_maf_file=${snp_maf_prefix}.refmtted.maf.tsv.gz

merged_alleles_file=${merged_alleles_dir}/01.merged.human_NHPs_alleles.tsv.gz

script=${SCRIPT_DIR}/merge_refmt_multizalign_with_chr_pos_allele.R
Rscript "$script" --snp_chr_pos_file "$human_alleles_file" --refmtted_multizalign_file "$refmtted_maf_file" \
    --out_file "$merged_alleles_file" 

### 4. find snps with human specific alleles
merged_alleles_with_HSAs=${merged_alleles_dir}/02.merged.human_NHPs_alleles.HSAs_added.tsv.gz

script=${SCRIPT_DIR}/find_human_specific_allele.py
python "$script" --input_file "$merged_alleles_file" --species_assembly_file "$species_assembly_file" \
    --out_file "$merged_alleles_with_HSAs"

### 5 Get SNPs with human specific allele and their position
## get snps with hsa
snps_with_HSAs=${merged_alleles_dir}/03.SNPs_with_HSAs.tsv.gz

gzip -dc "$merged_alleles_with_HSAs" | awk -F'\t' -vOFS='\t' '$58!="NA"{print $1,$2,$3,$4,$5,$58}' | bgzip > "$snps_with_HSAs"
tabix -f -S1 -s1 -b2 -e3 "$snps_with_HSAs"

## get snp pos
snp_pos_dir=${merged_alleles_dir}/snp_pos; mkdir -p "${snp_pos_dir}"

snps_with_HSAs_pos_with_chr_tsv=${snp_pos_dir}/snps_with_HSAs.pos.hg38.with_chr.tsv
snps_with_HSAs_pos_without_chr_tsv=${snp_pos_dir}/snps_with_HSAs.pos.hg38.without_chr.tsv

bgzip -dc "$snps_with_HSAs" | sed '1d' |cut -f1-3 > "$snps_with_HSAs_pos_with_chr_tsv"
sed 's/chr//' "$snps_with_HSAs_pos_with_chr_tsv" > "$snps_with_HSAs_pos_without_chr_tsv"

loginfo "Done. merged_alleles_dir is $merged_alleles_dir"
loginfo "SNPs with HSAs(Human Specific Alleles) is $snps_with_HSAs"

########################## Done ##########################
# 计时结束
endtime=$(date +"%Y-%m-%d %H:%M:%S")
"${SCRIPT_DIR}"/cal_used_time.sh "$starttime" "$endtime"
loginfo "Program is done"
