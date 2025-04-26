#!/usr/bin/env python
#%%
import os
import sys
import gzip
import time
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sparse

from ld_matrices_utils import load_alkes_ld, concat_ld_matrices, get_ld_for_snps

logging.basicConfig(level = logging.INFO, 
                    format = '%(asctime)s %(levelname)s %(filename)s %(funcName)s[line:%(lineno)d] : %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S'
                    )

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def run_shell_cmd(cmd):
    with subprocess.Popen(
        cmd, 
        shell=True, 
        bufsize=1, 
        executable = "/usr/bin/bash", 
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, 
        text=True
    ) as proc:
        for line in proc.stdout:
            print(line, end="")

    if proc.returncode != 0:
        sys.exit(proc.returncode)

def check_files(*args):
    for afile in args:
        if not os.path.exists(afile):
            raise Exception("Error! File `%s` not exists! Please check it!"%(afile))

def open_file(afile):
    if afile.endswith("gz"):
        file_reader = gzip.open(afile, "rt")
    else:
        file_reader = open(afile)
    return file_reader

def get_nrow(afile):
    count = 0
    file_reader = open_file(afile)
    for index, line in enumerate(file_reader):
        count += 1
    return count

def read_alkesgroup_smmr(smmr_file, header=True):
    '''
        header: SNP     CHR     POS     A1      A2      REF     EAF     Beta    se      P       N       INFO
        return: a iterator, with fields of 
            [SNP_ID, CHROM, POS, Effect_Allele, Other_Allele, Effect_Allele_Freq, BETA, SE, P, N, INFO]
    '''
    # if smmr_file.endswith('gz'):
    file_reader = open_file(smmr_file)
    if header:
        file_reader.readline()
    
    for l in file_reader.readlines():
        SNP, CHR, POS, A1, A2, REF, EAF, Beta, se, P, N, INFO = l.strip().split()
        yield SNP, CHR, POS, A1, A2, EAF, Beta, se, P, N, INFO

def get_ld_from_alkes(snp_to_extract_file, ld_prefixes, out_ld_file):
    '''
        从AlkesGroup提供的UKBB LD中获取与gwas summary相对应的LD信息
        snp_to_extract_file: Must contain at least four columns: chromosome, position, allele1 and allele2. Separated by tab. 
        ld_prefixes: A list of LD prefixes. 
        out_ld_file: Output file for extracted LD.

        Return: 
            A list contains SNPs which have LD. Format is chr_pos_ref_alt
    '''
    ## get snps (locus_allele) to find ld
    logging.info('Reading snp info from %s'%(snp_to_extract_file))
    snp_info = pd.read_csv(snp_to_extract_file, sep='\t')
    snp_info['locus_allele'] = snp_info[['chromosome', 'position', 'allele1', 'allele2']].astype(str).agg('_'.join, axis=1)
    snps_to_slct = snp_info['locus_allele'].tolist()

    ## load LD matrix
    logging.info('Loading LD matrix')
    if len(ld_prefixes) > 1:
        logging.info('There are multiple LD regions. Concatenating them as one...')
        ld_snps, ld_arr = concat_ld_matrices(ld_prefixes)

    else:
        first_ld_prefix = ld_prefixes[0]
        ld_snps, ld_arr = load_alkes_ld(first_ld_prefix)

    logging.info('Extracting LD matrix for SNPs')
    locus_alleles_with_ld, ld_for_snps = get_ld_for_snps(snps_to_slct, ld_snps, ld_arr)

    np.savetxt(out_ld_file, ld_for_snps, fmt='%.6f')
    logging.info(f'Saved LD into {out_ld_file}')

    return locus_alleles_with_ld

def extract_sumstats_for_locus(gwas_smmr_file, chrom, locus_start, locus_stop, out_sumstats_file):
    cmd = f"(bgzip -dc {gwas_smmr_file} | head -1; tabix {gwas_smmr_file} {chrom}:{locus_start}-{locus_stop}) | gzip > {out_sumstats_file}"
    logging.info('Getting sumstats for this locus, cmd is: \n%s'%(cmd))
    run_shell_cmd(cmd)

def get_sumstats_reader(sumstats_file, n_samples=None):
    # load sumstats
    smmr_iter = read_alkesgroup_smmr(sumstats_file, header=True)

    return smmr_iter

def filter_sumstats(sumstats_file, locus_start, locus_stop, 
                    min_info, min_freq, max_pval, 
                    out_maf_filtered_z_file, out_filtered_snp_file):
    '''
        Filter sumstats by MAF, INFO and P value
    '''
    out_maf_filtered_z_file_writer = open(out_maf_filtered_z_file, 'w')
    out_filtered_snp_file_writer = open(out_filtered_snp_file, 'w')
    num_filtered_snp = 0

    out_maf_filtered_z_file_writer.write('rsid\tchromosome\tposition\tallele1\tallele2\tmaf\tbeta\tse\n')
    out_filtered_snp_file_writer.write('SNP_ID\treason\n')

    sumstats_reader = get_sumstats_reader(sumstats_file)

    snplist = []
    for smmr in sumstats_reader:
        snp_id, chrom, pos, effect_allele, other_allele, effect_fllele_freq, beta, se, P, N, INFO = smmr
        if int(pos) >= locus_start and int(pos) <= locus_stop:

            # filter
            if float(INFO) < min_info:
                out_filtered_snp_file_writer.write('%s\tLow_INFO\n'%(snp_id))
                num_filtered_snp += 1
                continue

            # get freq
            freq = float(effect_fllele_freq)
            maf = freq if freq < 0.5 else 1 - freq
            maf = round(maf, 6)

            if maf < min_freq:
                out_filtered_snp_file_writer.write('%s\tMAF_too_low(%s)\n'%(snp_id, maf))
                out_filtered_snp_file_writer.flush()
                num_filtered_snp += 1
                continue
            
            if float(P) > max_pval:
                out_filtered_snp_file_writer.write('%s\tP_too_large(%s)\n'%(snp_id, P))
                out_filtered_snp_file_writer.flush()
                num_filtered_snp += 1
                continue

            # write results into z file
            out_maf_filtered_z_file_writer.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(
                snp_id, chrom, pos, effect_allele, other_allele, maf, beta, se
            ))
            out_maf_filtered_z_file_writer.flush()
            snplist.append(snp_id)

    out_maf_filtered_z_file_writer.close()
    out_filtered_snp_file_writer.close()
    logging.info('Removed %s SNPs due to low MAF, low INFO or large P value'%(num_filtered_snp))

    return snplist, num_filtered_snp

#%%
### load alkes LD
def prepare_fm_inputs(args):
    ## args
    chrom = args.chrom
    locus_start = args.locus_start
    locus_stop = args.locus_stop
    gwas_smmr_file = args.gwas_smmr
    ld_prefixes_file = args.ld_prefixes_file
    min_freq = args.min_freq
    min_info = args.min_info
    max_pval = args.max_pval
    out_dir = args.out_dir

    ## print args
    args_str = f'Arguments:\nchrom: {chrom}\nlocus_start: {locus_start}\nlocus_stop: {locus_stop}\n'
    args_str += f'gwas_smmr_file: {gwas_smmr_file}\n'
    args_str += f'ld_prefixes_file: {ld_prefixes_file}\nmin_freq: {min_freq}\nmin_info: {min_info}\n'
    args_str += f'max_pval: {max_pval}\nout_dir: {out_dir}\n\n'
    logging.info(args_str)

    ## Def output vars
    os.makedirs(out_dir, exist_ok=True)

    out_sumstats_dir = f'{out_dir}/01.sumstats'
    os.makedirs(out_sumstats_dir, exist_ok=True)

    filtered_snps_dir = f'{out_dir}/02.filtered_snps'
    os.makedirs(filtered_snps_dir, exist_ok=True)

    z_ld_dir = f'{out_dir}/03.z_ld'
    os.makedirs(z_ld_dir, exist_ok=True)

    ## Extract sumstats in this locus
    logging.info('1 / 4 | Extracting sumstats')
    out_sumstats_file = f'{out_sumstats_dir}/sumstats.tsv.gz'
    extract_sumstats_for_locus(gwas_smmr_file, chrom, locus_start, locus_stop, out_sumstats_file)

    ## Filter SNPs based on MAF, INFO and P value
    logging.info('2 / 4 | Filtering sumstats based on MAF, INFO and P value')
    out_maf_filtered_z_file = f'{filtered_snps_dir}/MAF_INFO_P_satisfied.z'
    out_filtered_snp_file = f'{filtered_snps_dir}/filtered_out_SNPs.txt'

    snplist, num_filtered_snp = filter_sumstats(out_sumstats_file, locus_start, locus_stop, 
                        min_info, min_freq, max_pval, 
                        out_maf_filtered_z_file, out_filtered_snp_file)

    ## Read LD prefixes
    ld_prefixes_file_reader = open_file(ld_prefixes_file)
    ld_prefixes = [x.strip() for x in ld_prefixes_file_reader]

    ## Get LD from UKBB LD provided by AlkesGroup
    out_ld_file = f'{z_ld_dir}/sumstats.ld'
    logging.info('3 / 4 | Extracting LD for SNPs')
    locus_alleles_with_ld = get_ld_from_alkes(out_maf_filtered_z_file, ld_prefixes, out_ld_file)

    ## Write sumstats of SNPs with LD to zfile
    logging.info('4 / 4 | Writing sumstats of SNPs with LD')
    out_z_file = f'{z_ld_dir}/sumstats.z'
    out_maf_filtered_z_file_reader = open_file(out_maf_filtered_z_file)
    out_z_file_writer = open(out_z_file, 'w')
    out_z_file_writer.write(out_maf_filtered_z_file_reader.readline())

    for l in out_maf_filtered_z_file_reader:
        line = l.strip().split()
        chrom, pos, allele1, allele2 = line[1:5]
        locus_allele = f'{chrom}_{pos}_{allele1}_{allele2}'

        if locus_allele in locus_alleles_with_ld:
            out_z_file_writer.write('%s\n'%('\t'.join(line)))

    ## Report summary
    num_total_snps = get_nrow(out_sumstats_file) - 1
    num_snps_for_finemapping = len(locus_alleles_with_ld)
    num_snps_without_ld = len(snplist) - num_snps_for_finemapping

    summary_str = f'Fine-mapping inputs are saved into {out_dir}.\nSummary of SNPs: \n'
    summary_str += f'Total number of SNPs in this locus: {num_total_snps}\n'
    summary_str += f'The number of SNPs removed due to low MAF, low INFO or large P value: {num_filtered_snp}\n'
    summary_str += f'The number of SNPs removed due to without LD: {num_snps_without_ld}\n'
    summary_str += f'The number of SNPs left for finemapping: {num_snps_for_finemapping}\n'

    logging.info(summary_str)

def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)

    parser.add_argument("--chrom", type=str, required=True, metavar='FILE', help="Chrom of lead SNP. No chr")
    parser.add_argument("--locus-start", type=int, required=True, metavar='FILE', help="Start pos of locus. ")
    parser.add_argument("--locus-stop", type=int, required=True, metavar='FILE', help="Stop pos of locus.")
    
    parser.add_argument("--gwas-smmr", type=str, required=True, metavar='FILE', help="GWAS summary statistics files. Must be compressed by bgzip and indexed by tabix.")
    parser.add_argument("--ld-prefixes-file", type=str, required=True, metavar='STR', help="File storing LD prefixes of UKBB LD from AlkesGroup. One line per LD prefix")
    parser.add_argument("--min-freq", type=float, default=0.0, metavar='FLOAT', help="Minimum allele frequency for SNPs to include [default 0]")
    parser.add_argument("--min-info", type=float, default=0.0, metavar='FLOAT', help="Minimum INFO score for SNPs to include [default 0]. Only used for AlkesGroup summary file, since summary of other two sources don't have INFO field")
    parser.add_argument("--max-pval", type=float, default=1, metavar='FLOAT', help="P value threshold [default 1]. Only SNPs with P less than this threshold can be used")
    parser.add_argument("--out-dir", type=str, required=True, metavar='STR', help="Output dir for input files.")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

def main():
    args = get_args()
    prepare_fm_inputs(args)

if __name__ == "__main__":
    main()
