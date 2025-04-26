#!/usr/bin/bash 

get_dirs(){
    # def inputs dirs
    inputs_dir=00.inputs; mkdir -p ${inputs_dir}
    gwas_sumstats_dir=${inputs_dir}/01.gwas_sumstats; mkdir -p ${gwas_sumstats_dir}
    ld_score_dir=${inputs_dir}/02.1KGP_Phase3_EUR_LDscore; mkdir -p ${ld_score_dir}
    lava_genomic_loci_dir=${inputs_dir}/03.lava_genomic_loci; mkdir -p ${lava_genomic_loci_dir}
    KGP_EUR_bfile_dir=${inputs_dir}/04.1KGP_EUR_bfile/; mkdir -p ${KGP_EUR_bfile_dir}
    UKBB_LD_dir=${inputs_dir}/05.UKBB_LD; mkdir -p ${UKBB_LD_dir}

    echo "Defining the following variables:

    inputs_dir=$inputs_dir
    gwas_sumstats_dir=$gwas_sumstats_dir
    ld_score_dir=$ld_score_dir
    lava_genomic_loci_dir=$lava_genomic_loci_dir
    KGP_EUR_bfile_dir=$KGP_EUR_bfile_dir
    UKBB_LD_dir=$UKBB_LD_dir
"
}


##################################################
########## 01. Calcluate global genetic correlation by LDSC
##################################################
ll scripts/wrkf/01.global_rg.sh


##################################################
########## 02. Calcluate local genetic correlation by LAVA
##################################################
ll scripts/wrkf/02.local_rg.sh


##################################################
########## 03. Perform statistical fine-mapping by SuSiE
##################################################
ll scripts/wrkf/03.finemapping.sh


##################################################
########## 04. Identify human-specific alleles (HSA) for the shared fine-mapped variants
##################################################
ll scripts/wrkf/04.hsa_for_fm_snps.sh


##################################################
########## 05. Perform functional analyses for shared fine-mapped variants with HSA
##################################################
ll scripts/wrkf/05.analyses_of_fm_snps_with_hsa.sh


