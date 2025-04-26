#!/usr/bin/bash 

get_dirs(){
    hsa_for_fm_snps_dir=04.hsa_for_fm_snps; mkdir -p ${hsa_for_fm_snps_dir}
    fm_snps_with_high_prob_hsa_dir=${hsa_for_fm_snps_dir}/04.fm_snps_with_high_prob_hsa; mkdir -p ${fm_snps_with_high_prob_hsa_dir}
    isec_fm_snps_pip_sum_gt05=${fm_snps_with_high_prob_hsa_dir}/03.HT_BMR_isec_FM_SNPS.with_HSAs.pip_sum_gt05.tsv.gz

    analyses_of_fm_snps_with_hsa_dir=05.analyses_of_fm_snps_with_hsa; mkdir -p ${analyses_of_fm_snps_with_hsa_dir}

    echo "Defining the following variables:
    hsa_for_fm_snps_dir=$hsa_for_fm_snps_dir
    fm_snps_with_high_prob_hsa_dir=$fm_snps_with_high_prob_hsa_dir
    isec_fm_snps_pip_sum_gt05=$isec_fm_snps_pip_sum_gt05

    # base dir
    analyses_of_fm_snps_with_hsa_dir=$analyses_of_fm_snps_with_hsa_dir
"
}

wrkf_script=scripts/wrkf/05.analyses_of_fm_snps_with_hsa.sh
# shellcheck source=/dev/null
source <(get-func.py -p $wrkf_script)
get_dirs


##################################################
########## 01. Prepare avinput format (01.avinput)
##################################################
### inputs
get_dirs

### outputs
avinput_dir=${analyses_of_fm_snps_with_hsa_dir}/01.avinput; mkdir -p ${avinput_dir}

### cmds
isec_fm_snps_pip_sum_gt05_avinput=${avinput_dir}/HT_BMR_isec_FM_SNPS.with_HSAs.pip_sum_gt05.avinput

bgzip -dc $isec_fm_snps_pip_sum_gt05 | sed '1d' | \
    awk -vOFS='\t' '{print "chr"$5, $6, $6, $7, $8}' > $isec_fm_snps_pip_sum_gt05_avinput


##################################################
########## 02. Annotate genes of these fm snps by ANNOVAR (02.annovar_gene_anno)
##################################################
### inputs
get_dirs

avinput_dir=${analyses_of_fm_snps_with_hsa_dir}/01.avinput; mkdir -p ${avinput_dir}
isec_fm_snps_pip_sum_gt05_avinput=${avinput_dir}/HT_BMR_isec_FM_SNPS.with_HSAs.pip_sum_gt05.avinput

### outputs
annovar_gene_anno_dir=${analyses_of_fm_snps_with_hsa_dir}/02.annovar_gene_anno; mkdir -p ${annovar_gene_anno_dir}
log_dir=${annovar_gene_anno_dir}/logs; mkdir -p ${log_dir}

### cmds
avinput=$isec_fm_snps_pip_sum_gt05_avinput
genome_version=hg19
threads=1
out_prefix=$annovar_gene_anno_dir/HT_BMR_isec_FM_SNPS
log=${log_dir}/perform_annovar_gene_anno.log

table_annovar="/path/to/annovar/table_annovar.pl"
humandb=/path/to/annovar/humandb/
$table_annovar -out "${out_prefix}" -buildver "${genome_version}" "${avinput}" "$humandb" \
    --protocol refGene -operation g -nastring . --thread "$threads" > $log 2>&1 &

# compress
multianno_file=${out_prefix}.${genome_version}_multianno.txt
bgzip "${multianno_file}"
tabix -S1 -s1 -b2 -e2 "${multianno_file}".gz

# get region counts
multianno_file=$annovar_gene_anno_dir/HT_BMR_isec_FM_SNPS.${genome_version}_multianno.txt.gz
region_counts=${annovar_gene_anno_dir}/region_counts.tsv

(
    echo -e "region\tnum_snps";
    bgzip -dc $multianno_file |sed '1d' |cut -f6 |sort |uniq -c |awk -vOFS='\t' '{print $2,$1}'
) > $region_counts

cat $region_counts
# region  num_snps
# downstream      30
# exonic  82
# intergenic      1052
# intronic        1715
# ncRNA_exonic    17
# ncRNA_intronic  277
# splicing        1
# upstream        38
# upstream;downstream     3
# UTR3    79
# UTR5    12
# UTR5;UTR3       1


