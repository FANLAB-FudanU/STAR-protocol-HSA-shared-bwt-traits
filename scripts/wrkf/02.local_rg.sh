#!/usr/bin/bash 

########## Main: 02. Calcluate local genetic correlation by LAVA

get_dirs(){
    # def inputs dirs
    inputs_dir=00.inputs; mkdir -p ${inputs_dir}
    gwas_sumstats_dir=${inputs_dir}/01.gwas_sumstats; mkdir -p ${gwas_sumstats_dir}
    ld_score_dir=${inputs_dir}/02.1KGP_Phase3_EUR_LDscore; mkdir -p ${ld_score_dir}
    lava_genomic_loci_dir=${inputs_dir}/03.lava_genomic_loci; mkdir -p ${lava_genomic_loci_dir}
    KGP_EUR_bfile_dir=${inputs_dir}/04.1KGP_EUR_bfile/; mkdir -p ${KGP_EUR_bfile_dir}

    local_rg_dir=02.local_rg; mkdir -p ${local_rg_dir}
    
    echo "Defining the following variables:

    inputs_dir=$inputs_dir
    gwas_sumstats_dir=$gwas_sumstats_dir
    ld_score_dir=$ld_score_dir
    lava_genomic_loci_dir=$lava_genomic_loci_dir
    KGP_EUR_bfile_dir=$KGP_EUR_bfile_dir

    # base dir
    local_rg_dir=$local_rg_dir
"
}


##################################################
########## 01. Prepare inputs for LAVA (01.lava_inputs)
##################################################
### inputs
get_dirs

### outputs
lava_inputs_dir=${local_rg_dir}/01.lava_inputs; mkdir -p ${lava_inputs_dir}
log_dir=${lava_inputs_dir}/logs; mkdir -p ${log_dir}

### cmds
trait_ids=(HT BMR)
min_info=0.9
min_maf=0.01

script=scripts/bin/make_lava_format_sumstats.sh
for trait_id in "${trait_ids[@]}"
do
    gwas_file=${gwas_sumstats_dir}/${trait_id}.sumstats.tsv.gz
    out_file=${lava_inputs_dir}/${trait_id}.LAVA_fmt.tsv.gz
    log=${log_dir}/prepare_inputs_for_${trait_id}.log

    bash $script "$gwas_file" "$min_info" $min_maf "$out_file" > "$log" 2>&1
    
    echo "done for $trait_id"
done


##################################################
########## 02. Make sample overlap file from LDSC outputs (02.sample_overlap)
##################################################
## 这一步是参考了LAVA给的代码 LAVA/vignettes/sample_overlap.md
### inputs
get_dirs

global_rg_dir=01.global_rg; mkdir -p ${global_rg_dir}
rg_dir=${global_rg_dir}/02.rg; mkdir -p ${rg_dir}
HT_BMR_rg_summary=${rg_dir}/rg_summary.tsv

### outputs
sample_overlap_dir=${local_rg_dir}/02.sample_overlap; mkdir -p ${sample_overlap_dir}
log_dir=${local_rg_dir}/logs; mkdir -p ${log_dir}

### cmds
sample_overlap_file=${sample_overlap_dir}/HT_BMR.sample_overlap.tsv
log=${log_dir}/create_HT_BMR_sample_overlap.log

export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
script=scripts/bin/create_sample_overlap.R
Rscript $script --rg_summary $HT_BMR_rg_summary --out_file $sample_overlap_file > $log 2>&1 &


##################################################
########## 03. Run LAVA using whole-genome wide SNPs in 2495 loci (03.lava_outputs)
##################################################
### inputs
get_dirs

lava_inputs_dir=${local_rg_dir}/01.lava_inputs; mkdir -p ${lava_inputs_dir}
HT_LAVA_sumstats=${lava_inputs_dir}/HT.LAVA_fmt.tsv.gz
BMR_LAVA_sumstats=${lava_inputs_dir}/BMR.LAVA_fmt.tsv.gz

lava_loc_file=${lava_genomic_loci_dir}/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile

eur_bfile=${KGP_EUR_bfile_dir}/EUR

sample_overlap_dir=${local_rg_dir}/02.sample_overlap; mkdir -p ${sample_overlap_dir}
sample_overlap_file=${sample_overlap_dir}/HT_BMR.sample_overlap.tsv

### outputs
lava_outputs_dir=${local_rg_dir}/03.lava_outputs; mkdir -p ${lava_outputs_dir}
log_dir=${local_rg_dir}/logs; mkdir -p ${log_dir}

### cmds
## Generate input info file
input_info_file=${lava_outputs_dir}/01.HT_BMR_input_info_file.tsv

(
    echo -e "phenotype\tcases\tcontrols\tfilename"; 
    echo -e "HT\tNA\tNA\t$HT_LAVA_sumstats";
    echo -e "BMR\tNA\tNA\t$BMR_LAVA_sumstats"
) > $input_info_file

## run lava
input_info_file=${lava_outputs_dir}/01.HT_BMR_input_info_file.tsv
ref_bfile=$eur_bfile
loci_file=$lava_loc_file
threads=40
out_prefix=${lava_outputs_dir}/02.HT_BMR_local_rg
log=$log_dir/run_lava_HT_BMR.log

export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
script=scripts/bin/run_lava.R
Rscript $script --input_info_file "$input_info_file" --sample_overlap_file "$sample_overlap_file" --ref_bfile "$ref_bfile" \
    --loci_file "$loci_file" --threads $threads --out_prefix "$out_prefix" > $log 2>&1 &

## plot local rg
sig_local_rg_file=${lava_outputs_dir}/02.HT_BMR_local_rg.sig_bivar.tsv
title='Distribution of local rg across the genome'
unit_mm=TRUE
width=190
height=180
rg_heatmap_file=${lava_outputs_dir}/03.local_rg_heatmap.pdf
log=${log_dir}/plot_genomic_heatmap_of_local_rg.log

script=scripts/bin/plot_local_rg.R
export PATH=/path/to/miniforge3/envs/r4.3/bin/:$PATH
Rscript $script --sig_local_rg_file $sig_local_rg_file \
    --title "$title" --unit_mm $unit_mm --width $width --height $height \
    --out_file $rg_heatmap_file > $log 2>&1 &



