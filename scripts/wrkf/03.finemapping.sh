#!/usr/bin/bash 

get_dirs(){
    # def inputs dirs
    inputs_dir=00.inputs; mkdir -p ${inputs_dir}
    gwas_sumstats_dir=${inputs_dir}/01.gwas_sumstats; mkdir -p ${gwas_sumstats_dir}
    ld_score_dir=${inputs_dir}/02.1KGP_Phase3_EUR_LDscore; mkdir -p ${ld_score_dir}
    lava_genomic_loci_dir=${inputs_dir}/03.lava_genomic_loci; mkdir -p ${lava_genomic_loci_dir}
    KGP_EUR_bfile_dir=${inputs_dir}/04.1KGP_EUR_bfile/; mkdir -p ${KGP_EUR_bfile_dir}
    UKBB_LD_dir=${inputs_dir}/05.UKBB_LD; mkdir -p ${UKBB_LD_dir}

    finemapping_dir=03.finemapping; mkdir -p ${finemapping_dir}

    echo "Defining the following variables:

    inputs_dir=$inputs_dir
    gwas_sumstats_dir=$gwas_sumstats_dir
    ld_score_dir=$ld_score_dir
    lava_genomic_loci_dir=$lava_genomic_loci_dir
    KGP_EUR_bfile_dir=$KGP_EUR_bfile_dir
    UKBB_LD_dir=$UKBB_LD_dir

    # base dir
    finemapping_dir=$finemapping_dir
"
}

##################################################
########## 01. Get loci with significant local rg and add ld_prefix to loci file (01.loci_with_sig_rg)
##################################################
### inputs
get_dirs

local_rg_dir=02.local_rg; mkdir -p ${local_rg_dir}
lava_outputs_dir=${local_rg_dir}/03.lava_outputs; mkdir -p ${lava_outputs_dir}
lava_sig_bivar_file=${lava_outputs_dir}/02.HT_BMR_local_rg.sig_bivar.tsv

### outputs
loci_with_sig_rg_dir=${finemapping_dir}/01.loci_with_sig_rg; mkdir -p ${loci_with_sig_rg_dir}
log_dir=${loci_with_sig_rg_dir}/logs; mkdir -p ${log_dir}

### cmds
## get loci with sig local rg
loci_with_sig_local_rg=${loci_with_sig_rg_dir}/01.loci_with_sig_local_rg.tsv

awk -vOFS='\t' '{print $1,$2,$3,$4,$5,$9,$12,$15}' $lava_sig_bivar_file > ${loci_with_sig_local_rg}

## add ld prefix to loci with sig local rg
ld_dir=${UKBB_LD_dir}
loci_file=${loci_with_sig_local_rg}
out_prefix=${loci_with_sig_rg_dir}/02.loci_with_sig_local_rg
log=${log_dir}/find_ld_prefix_of_loci_with_sig_local_rg.log

export PATH=/path/to/miniforge3/envs/py3.10/bin/:$PATH
script=scripts/bin/get_ld_prefix_for_loci.sh
bash $script "$ld_dir" "$loci_file" "$out_prefix" > $log 2>&1 &

## we only perform fine-mapping for loci in ${loci_with_sig_rg_dir}/02.loci_with_sig_local_rg.with_npz_ld.tsv
# since these loci do not contain long-range LD


##################################################
########## 02. Perform fine-mapping for height (HT) and BMR, respectively (02.per_locus_fm)
##################################################
### inputs
get_dirs

work_dir=/path/to/HSA-protocol

loci_with_sig_rg_dir=${finemapping_dir}/01.loci_with_sig_rg; mkdir -p ${loci_with_sig_rg_dir}
loci_with_sig_local_rg_ld_prefix=${loci_with_sig_rg_dir}/02.loci_with_sig_local_rg.with_npz_ld.tsv

### outputs
per_locus_fm_dir=${finemapping_dir}/02.per_locus_fm; mkdir -p ${per_locus_fm_dir}
log_dir=${per_locus_fm_dir}/logs; mkdir -p ${log_dir}

### cmds
loci_file=${loci_with_sig_local_rg_ld_prefix}
min_freq=0.01
min_info=0.9
max_pval=1
n_samples=458303
n_causal_snps=10
threads=20
dry_run=False

trait_ids=(HT BMR)

export PATH=/path/to/miniforge3/envs/py3.10/bin/:$PATH
export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
script=scripts/bin/perform_fm_for_loci_file.sh
for trait_id in "${trait_ids[@]}"
do
    gwas_sumstats=$gwas_sumstats_dir/${trait_id}.sumstats.tsv.gz
    out_dir=${per_locus_fm_dir}/${trait_id}
    log=${log_dir}/perform_fm_with_${trait_id}_sumstats.log

    bash $script "$work_dir" "$loci_file" "$gwas_sumstats" "$min_freq" "$min_info" \
        "$max_pval" "$n_samples" "$n_causal_snps" "$threads" "$dry_run" "$out_dir" > "$log" 2>&1
    echo "done for $trait_id"
done &

## timing: ~1 hour for each trait using 20 threads


##################################################
########## 03. Get stats of fine-mapped SNPs (fm-snps) for HT and BMR (03.fm_snps_stats)
##################################################
### inputs
get_dirs

per_locus_fm_dir=${finemapping_dir}/02.per_locus_fm; mkdir -p ${per_locus_fm_dir}

loci_with_sig_rg_dir=${finemapping_dir}/01.loci_with_sig_rg; mkdir -p ${loci_with_sig_rg_dir}
loci_with_sig_local_rg_ld_prefix=${loci_with_sig_rg_dir}/02.loci_with_sig_local_rg.with_npz_ld.tsv

### outputs
fm_snps_stats_dir=${finemapping_dir}/03.fm_snps_stats; mkdir -p ${fm_snps_stats_dir}
log_dir=${fm_snps_stats_dir}/logs; mkdir -p ${log_dir}

### cmds
loci_file=${loci_with_sig_local_rg_ld_prefix}

trait_ids=(HT BMR)
export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
script=scripts/bin/summarize_fm_for_loci_file.sh
for trait_id in "${trait_ids[@]}"
do
    trait_per_locus_fm_dir=${per_locus_fm_dir}/${trait_id}/per_locus
    out_dir=${fm_snps_stats_dir}/${trait_id}
    log=${log_dir}/get_${trait_id}_fm_stats.log

    bash $script "$trait_per_locus_fm_dir" "$loci_file" "$out_dir" > "$log" 2>&1
    echo "done for $trait_id"
done


##################################################
########## 04. Merge fm snps of all loci (04.all_loci_merged_fm_snps)
##################################################
### inputs
get_dirs

per_locus_fm_dir=${finemapping_dir}/02.per_locus_fm; mkdir -p ${per_locus_fm_dir}
fm_snps_stats_dir=${finemapping_dir}/03.fm_snps_stats; mkdir -p ${fm_snps_stats_dir}

### outputs
all_loci_merged_fm_snps_dir=${finemapping_dir}/04.all_loci_merged_fm_snps; mkdir -p ${all_loci_merged_fm_snps_dir}
log_dir=${all_loci_merged_fm_snps_dir}/logs; mkdir -p ${log_dir}

### cmds
trait_ids=(HT BMR)

export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
script=scripts/bin/merge_per_locus_fm_snps.sh
for trait_id in "${trait_ids[@]}"
do
    trait_per_locus_fm_dir=${per_locus_fm_dir}/${trait_id}/per_locus
    loci_with_fm_snps=${fm_snps_stats_dir}/${trait_id}/loci_with_fm_snps_stats.tsv
    out_dir=${all_loci_merged_fm_snps_dir}/${trait_id}
    log=${log_dir}/merge_${trait_id}_fm_snps.log

    bash $script "$trait_per_locus_fm_dir" "$loci_with_fm_snps" "$out_dir" > "$log" 2>&1 &
    echo "done for $trait_id"
done


##################################################
########## 05. Get intersected fm snps of HT and BMR (05.HT_BMR_isec_fm_snps)
##################################################
### inputs
get_dirs

all_loci_merged_fm_snps_dir=${finemapping_dir}/04.all_loci_merged_fm_snps; mkdir -p ${all_loci_merged_fm_snps_dir}
HT_merged_fm_snps=${all_loci_merged_fm_snps_dir}/HT/merged_FM_SNPs.tsv.gz
BMR_merged_fm_snps=${all_loci_merged_fm_snps_dir}/BMR/merged_FM_SNPs.tsv.gz

### outputs
HT_BMR_isec_fm_snps_dir=${finemapping_dir}/05.HT_BMR_isec_fm_snps; mkdir -p ${HT_BMR_isec_fm_snps_dir}

### cmds
isec_fm_snps_file=${HT_BMR_isec_fm_snps_dir}/HT_BMR_isec_FM_SNPs.tsv.gz

csvtk join -t -T $HT_merged_fm_snps $BMR_merged_fm_snps -f "1,2,3,4,5,6,7,8" -s "HT,BMR" | \
    bgzip > $isec_fm_snps_file
tabix -S1 -s5 -b6 -e6 $isec_fm_snps_file



