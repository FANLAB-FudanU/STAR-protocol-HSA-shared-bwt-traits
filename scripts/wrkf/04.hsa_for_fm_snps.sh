#!/usr/bin/bash 

get_dirs(){
    # def inputs dirs
    inputs_dir=00.inputs; mkdir -p ${inputs_dir}
    gwas_sumstats_dir=${inputs_dir}/01.gwas_sumstats; mkdir -p ${gwas_sumstats_dir}
    ld_score_dir=${inputs_dir}/02.1KGP_Phase3_EUR_LDscore; mkdir -p ${ld_score_dir}
    lava_genomic_loci_dir=${inputs_dir}/03.lava_genomic_loci; mkdir -p ${lava_genomic_loci_dir}
    KGP_EUR_bfile_dir=${inputs_dir}/04.1KGP_EUR_bfile/; mkdir -p ${KGP_EUR_bfile_dir}
    UKBB_LD_dir=${inputs_dir}/05.UKBB_LD; mkdir -p ${UKBB_LD_dir}
    UCSC_dbSNP_bigbed_dir=${inputs_dir}/06.UCSC_dbSNP_bigbed/; mkdir -p ${UCSC_dbSNP_bigbed_dir}

    finemapping_dir=03.finemapping; mkdir -p ${finemapping_dir}

    hsa_for_fm_snps_dir=04.hsa_for_fm_snps; mkdir -p ${hsa_for_fm_snps_dir}

    echo "Defining the following variables:

    inputs_dir=$inputs_dir
    gwas_sumstats_dir=$gwas_sumstats_dir
    ld_score_dir=$ld_score_dir
    lava_genomic_loci_dir=$lava_genomic_loci_dir
    KGP_EUR_bfile_dir=$KGP_EUR_bfile_dir
    UKBB_LD_dir=$UKBB_LD_dir
    UCSC_dbSNP_bigbed_dir=$UCSC_dbSNP_bigbed_dir

    finemapping_dir=$finemapping_dir

    # base dir
    hsa_for_fm_snps_dir=$hsa_for_fm_snps_dir
"
}

##################################################
########## 01. Get hg38 of fm snps (01.fm_snps_in_hg38)
##################################################
### inputs
get_dirs

work_dir=/path/to/HSA-protocol
dbsnp_hg19_bigbed_file=${UCSC_dbSNP_bigbed_dir}/dbSnp155.hg19.bb
dbsnp_hg38_bigbed_file=${UCSC_dbSNP_bigbed_dir}/dbSnp155.hg38.bb

HT_BMR_isec_fm_snps_dir=${finemapping_dir}/05.HT_BMR_isec_fm_snps; mkdir -p ${HT_BMR_isec_fm_snps_dir}
HT_BMR_isec_fm_snps_file=${HT_BMR_isec_fm_snps_dir}/HT_BMR_isec_FM_SNPs.tsv.gz

### outputs
fm_snps_in_hg38_dir=${hsa_for_fm_snps_dir}/01.fm_snps_in_hg38; mkdir -p ${fm_snps_in_hg38_dir}
log_dir=${fm_snps_in_hg38_dir}/logs; mkdir -p ${log_dir}

### cmds
## get rsid of fm snps
fm_snps_rsid=${fm_snps_in_hg38_dir}/01.fm_snps.rsid.txt

bgzip -dc $HT_BMR_isec_fm_snps_file | sed '1d' | cut -f4 | sort -u > $fm_snps_rsid

## generate chroms list
chroms_list=${fm_snps_in_hg38_dir}/02.chroms_list.txt

bgzip -dc $HT_BMR_isec_fm_snps_file |sed '1d' |awk '{print "chr"$1}' |sort -u > $chroms_list

## get hg19 and hg38 position of fm snps
rsid_file=$fm_snps_rsid
chroms_list=${fm_snps_in_hg38_dir}/02.chroms_list.txt
num_rows=1000
threads=4
out_dir=$fm_snps_in_hg38_dir/03.fm_snps_in_hg38_pos
log=${log_dir}/get_hg19_hg38_pos_of_isec_fm_snps.log

script=scripts/bin/get_hg19_hg38_snp_pos_by_rsid.sh
export PATH=/path/to/miniforge3/envs/py3.10/bin/:$PATH
export PATH=/path/to/ucsc_utility/:$PATH
bash $script "$work_dir" "$dbsnp_hg38_bigbed_file" "$dbsnp_hg19_bigbed_file" "$rsid_file" \
    "$chroms_list" "$num_rows" "$threads" "$out_dir" > $log 2>&1 &


##################################################
########## 02. Get hg38 pos and alleles for fm snps (02.fm_snps_pos_alleles_in_hg38)
##################################################
### inputs
get_dirs

# HT BMR isec fm snps
HT_BMR_isec_fm_snps_dir=${finemapping_dir}/05.HT_BMR_isec_fm_snps; mkdir -p ${HT_BMR_isec_fm_snps_dir}
HT_BMR_isec_fm_snps_file=${HT_BMR_isec_fm_snps_dir}/HT_BMR_isec_FM_SNPs.tsv.gz

# hg19 to hg38 coord lookup table of HT BMR isec fm snps
fm_snps_in_hg38_dir=${hsa_for_fm_snps_dir}/01.fm_snps_in_hg38; mkdir -p ${fm_snps_in_hg38_dir}
fm_snps_hg19_to_hg38_coord_lookup_table=$fm_snps_in_hg38_dir/03.fm_snps_in_hg38_pos/03.coord_lookup_table/hg19_to_hg38_coord_lookup_table.with_chr.tsv.gz

### outputs
fm_snps_pos_alleles_in_hg38_dir=${hsa_for_fm_snps_dir}/02.fm_snps_pos_alleles_in_hg38; mkdir -p ${fm_snps_pos_alleles_in_hg38_dir}
log_dir=${hsa_for_fm_snps_dir}/logs; mkdir -p ${log_dir}

### cmds
hg19_to_hg38_coord_lookup_table=$fm_snps_hg19_to_hg38_coord_lookup_table
isec_fm_snps_file=$HT_BMR_isec_fm_snps_file
out_file=${fm_snps_pos_alleles_in_hg38_dir}/fm_snps_hg38_pos_alleles.tsv
log=${log_dir}/get_hg38_pos_and_alleles_for_fm_snps.log

script=scripts/bin/prepare_hsa_inputs.sh
bash $script --hg19_to_hg38_coord_lookup_table="$hg19_to_hg38_coord_lookup_table" \
    --isec_fm_snps_file="$isec_fm_snps_file" --out_file="$out_file" > $log 2>&1 &


##################################################
########## 03. Identify human-specific alleles for fm snps (03.human_specific_alleles)
##################################################
### inputs
get_dirs

work_dir=/path/to/HSA-protocol

fm_snps_pos_alleles_in_hg38_dir=${hsa_for_fm_snps_dir}/02.fm_snps_pos_alleles_in_hg38; mkdir -p ${fm_snps_pos_alleles_in_hg38_dir}
fm_snps_hg38_pos_alleles=${fm_snps_pos_alleles_in_hg38_dir}/fm_snps_hg38_pos_alleles.tsv

fm_snps_in_hg38_dir=${hsa_for_fm_snps_dir}/01.fm_snps_in_hg38; mkdir -p ${fm_snps_in_hg38_dir}
fm_snps_hg38_to_hg19_coord_lookup_table=$fm_snps_in_hg38_dir/03.fm_snps_in_hg38_pos/03.coord_lookup_table/hg38_to_hg19_coord_lookup_table.with_chr.tsv.gz

ucsc_multiz30way_dir=${inputs_dir}/06.UCSC_multiz30way
species_assembly_file=${inputs_dir}/07.cons30primates_hg38.non_human_primates.lookup_table.tsv

### outputs
human_specific_alleles_dir=${hsa_for_fm_snps_dir}/03.human_specific_alleles; mkdir -p ${human_specific_alleles_dir}
log_dir=${human_specific_alleles_dir}/logs; mkdir -p ${log_dir}

### cmds
## identify hsa for fm snps (hg38) (5 mins)
human_alleles_file=${fm_snps_hg38_pos_alleles}
raw_maf_dir=${ucsc_multiz30way_dir}
human_assembly=hg38
threads=22
out_dir=${human_specific_alleles_dir}/01.hsa_in_hg38
log=${log_dir}/get_isec_fm_snps_with_hsa.log

script=scripts/bin/get_snps_with_human_specific_alleles.sh
export PATH=/path/to/miniforge3/envs/py3.10/bin/:$PATH
bash $script "$work_dir" "$human_alleles_file" "$raw_maf_dir" "$human_assembly" "$species_assembly_file" \
    "$threads" "$out_dir" > $log 2>&1 &

## convert hsa for fm snps to hg19
isec_fm_snps_with_hsa_hg38_file=${human_specific_alleles_dir}/01.hsa_in_hg38/02.merged_alleles/03.SNPs_with_HSAs.tsv.gz
isec_fm_snps_with_hsa_hg19_file=${human_specific_alleles_dir}/02.isec_fm_snps_with_HSA.hg19.tsv
log=${log_dir}/convert_hsa_to_hg19.log

script=scripts/bin/convert_hsa_snps_genomic_versions.py
export PATH=/path/to/miniforge3/envs/py3.10/bin/:$PATH
python -u $script --input_file $isec_fm_snps_with_hsa_hg38_file --coord_lookup_table_file $fm_snps_hg38_to_hg19_coord_lookup_table \
    --out_file $isec_fm_snps_with_hsa_hg19_file > $log 2>&1 &

## get position of fm snps of hsa in hg19
isec_fm_snps_with_hsa_hg19_file=${human_specific_alleles_dir}/02.isec_fm_snps_with_HSA.hg19.tsv.gz

isec_fm_snps_with_hsa_hg19_pos_without_chr=${human_specific_alleles_dir}/03.isec_fm_snps_with_HSA.hg19.pos.without_chr.tsv

bgzip -dc $isec_fm_snps_with_hsa_hg19_file | sed '1d' | awk -vOFS='\t' '{print $1,$2,$3}' |sed 's/chr//' > $isec_fm_snps_with_hsa_hg19_pos_without_chr


##################################################
########## 04. Identify fm SNPs with HSA with high probability (04.fm_snps_with_high_prob_hsa)
##################################################
# Now we need to identify credible sets whose fm SNPs have a total of PIP >= 50%
# To do this, we need to extract fine-mapping results of fm SNPs with HSA.
# Then we calculate the total PIP of each credible sets. 
# Finally we keep credible sets whose total PIP are >= 50%
### inputs
get_dirs

# HT BMR isec fm snps
HT_BMR_isec_fm_snps_dir=${finemapping_dir}/05.HT_BMR_isec_fm_snps; mkdir -p ${HT_BMR_isec_fm_snps_dir}
HT_BMR_isec_fm_snps_file=${HT_BMR_isec_fm_snps_dir}/HT_BMR_isec_FM_SNPs.tsv.gz

# hg19 pos of isec fm snps with hsa
isec_fm_snps_with_hsa_hg19_pos_without_chr=${human_specific_alleles_dir}/03.isec_fm_snps_with_HSA.hg19.pos.without_chr.tsv

### outputs
fm_snps_with_high_prob_hsa_dir=${hsa_for_fm_snps_dir}/04.fm_snps_with_high_prob_hsa; mkdir -p ${fm_snps_with_high_prob_hsa_dir}
log_dir=${fm_snps_with_high_prob_hsa_dir}/logs; mkdir -p ${log_dir}

### cmds
## get fine-mapping results of isec fm snps with hsa
isec_fm_snps_with_hsa_fm_res_file=${fm_snps_with_high_prob_hsa_dir}/01.HT_BMR_isec_FM_SNPs.with_HSAs.tsv.gz
(
    bgzip -dc $HT_BMR_isec_fm_snps_file | head -1;
    tabix $HT_BMR_isec_fm_snps_file -R $isec_fm_snps_with_hsa_hg19_pos_without_chr
) | bgzip > $isec_fm_snps_with_hsa_fm_res_file
tabix -S1 -s5 -b6 -e6 $isec_fm_snps_with_hsa_fm_res_file

## calculate the total PIP of each credible sets
isec_fm_snps_with_hsa_fm_res_file=${fm_snps_with_high_prob_hsa_dir}/01.HT_BMR_isec_FM_SNPs.with_HSAs.tsv.gz
isec_fm_snps_with_hsa_pip_sum_dir=${fm_snps_with_high_prob_hsa_dir}/02.HT_BMR_isec_FM_SNPs.with_HSAs.sum_of_pips_per_cs
log=${log_dir}/get_sum_of_pips_of_isec_fm_snps.log

script=scripts/bin/cal_sum_of_pip_for_isec_fm_snps.R
export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
Rscript $script --isec_fm_snps_file $isec_fm_snps_with_hsa_fm_res_file \
    --out_dir $isec_fm_snps_with_hsa_pip_sum_dir > $log 2>&1 &

## keep credible sets whose total PIP are >= 50%
isec_fm_snps_with_hsa_pip_sum_dir=${fm_snps_with_high_prob_hsa_dir}/02.HT_BMR_isec_FM_SNPs.with_HSAs.sum_of_pips_per_cs
cs_with_pip_gt05=${isec_fm_snps_with_hsa_pip_sum_dir}/cs_with_pip_gt05.HT.tsv

isec_fm_snps_pip_sum_gt05_fm_res=${fm_snps_with_high_prob_hsa_dir}/03.HT_BMR_isec_FM_SNPS.with_HSAs.pip_sum_gt05.tsv.gz

bgzip -dc $isec_fm_snps_with_hsa_fm_res_file | awk -F'\t' -vOFS='\t' 'NR==FNR{
    cs_id=$1"_"$2"_"$3"_"$4;
    cs_ids[cs_id]=1;
} NR>FNR {
    cs_id=$1"_"$2"_"$3"_"$15;
    if(cs_id in cs_ids){
        print $0
    }
}' $cs_with_pip_gt05 - | bgzip > $isec_fm_snps_pip_sum_gt05_fm_res
tabix -S1 -s5 -b6 -e6 $isec_fm_snps_pip_sum_gt05_fm_res


