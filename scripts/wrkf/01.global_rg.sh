#!/usr/bin/bash 

########## Main: 01. Calcluate global genetic correlation by LDSC

get_dirs(){
    # def inputs dirs
    inputs_dir=00.inputs; mkdir -p ${inputs_dir}
    gwas_sumstats_dir=${inputs_dir}/01.gwas_sumstats; mkdir -p ${gwas_sumstats_dir}
    ld_score_dir=${inputs_dir}/02.1KGP_Phase3_EUR_LDscore; mkdir -p ${ld_score_dir}
    
    global_rg_dir=01.global_rg; mkdir -p ${global_rg_dir}
    
    echo "Defining the following variables:

    inputs_dir=$inputs_dir
    gwas_sumstats_dir=$gwas_sumstats_dir
    ld_score_dir=$ld_score_dir

    # base dir
    global_rg_dir=$global_rg_dir
"
}

##################################################
########## 01. Convert gwas sumstats to LDSC input files (01.ldsc_inputs)
##################################################
### inputs
get_dirs

### outputs
ldsc_inputs_dir=${global_rg_dir}/01.ldsc_inputs; mkdir -p ${ldsc_inputs_dir}
log_dir=${ldsc_inputs_dir}/logs; mkdir -p ${log_dir}

### cmds
trait_ids=(HT BMR)

export PATH=/path/to/miniforge3/envs/py2.7/bin/:$PATH
munge_sumstats_script=/path/to/ldsc/munge_sumstats.py
for trait_id in "${trait_ids[@]}"
do
    gwas_file=${gwas_sumstats_dir}/${trait_id}.sumstats.tsv.gz
    out_prefix=${ldsc_inputs_dir}/${trait_id}
    log=${log_dir}/prepare_inputs_for_${trait_id}.log

    python -u $munge_sumstats_script --sumstats "$gwas_file" \
        --out "${out_prefix}" > "$log"
    
    echo "done for $trait_id"
done


##################################################
########## 02. Calculate genetic correlation (02.rg)
##################################################
### inputs
get_dirs

ldsc_inputs_dir=${global_rg_dir}/01.ldsc_inputs; mkdir -p ${ldsc_inputs_dir}
ld_score_str=$ld_score_dir/LDscore.

### outputs
rg_dir=${global_rg_dir}/02.rg; mkdir -p ${rg_dir}
log_dir=${global_rg_dir}/logs; mkdir -p ${log_dir}

### cmds
## cal rg
ldsc_input_files=${ldsc_inputs_dir}/HT.sumstats.gz,${ldsc_inputs_dir}/BMR.sumstats.gz
out_prefix=${rg_dir}/HT_BMR
log=${log_dir}/cal_rg_bwt_HT_BMR.log

export PATH=/path/to/miniforge3/envs/py2.7/bin/:$PATH
ldsc_script=/path/to/ldsc/ldsc.py
python -u $ldsc_script \
    --rg "$ldsc_input_files" \
    --ref-ld-chr "$ld_score_str" \
    --w-ld-chr "$ld_score_str" \
    --out "$out_prefix" > $log 2>&1 &

## reformat output file
out_file=${rg_dir}/rg_summary.tsv
(
    grep -A 1 'Summary of Genetic Correlation Results' "$log" | sed '1d' \
        | awk '{for(i=1;i<NF;i++){printf("%s\t", $i)}; print $NF }';
    grep -A 2 'Summary of Genetic Correlation Results' "$log" | tail -n1 \
        | awk '{printf("HT\tBMR\t"); for(i=3;i<NF;i++){printf("%s\t", $i)}; print $NF }'
) > "$out_file"




