#!/usr/bin/env python
#%%
import os
import sys
import gzip
import time
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sparse

logging.basicConfig(level = logging.INFO, 
                    format = '%(asctime)s %(levelname)s %(filename)s %(funcName)s[line:%(lineno)d] : %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S'
                    )

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """
        可以将参数的默认选项打印在帮助文档中
    """
    pass

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

def load_alkes_ld(ld_prefix):
    '''
        ld_prefix = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr10_135000001_138000001'
        ld_snps, ld_arr = load_alkes_ld(ld_prefix)
    '''
    # get snp
    ld_file = '%s.npz' % (ld_prefix)
    ld_file2 = '%s.npz2' % (ld_prefix)
    ld_snp_file = '%s.gz' % (ld_prefix)
    if not os.path.exists(ld_snp_file):
        raise Exception(f'Error! {ld_snp_file} does not exist!')

    logging.info('Reading SNPs from: %s'%(ld_snp_file))
    ld_snps = pd.read_csv(ld_snp_file, sep='\t')

    # get LD
    if os.path.exists(ld_file):
        logging.info('Reading LD from npz file: %s'%(ld_file))
        ld_arr = sparse.load_npz(ld_file).toarray()

    elif os.path.exists(ld_file2):
        logging.info('Reading LD from npz2 file: %s'%(ld_file2))
        ld_arr = sparse.load_npz(ld_file2).toarray()

    else:
        raise Exception(f'Error! Both {ld_file} and {ld_file2} not exist!')

    ld_arr = ld_arr + ld_arr.T
    assert np.allclose(np.diag(ld_arr), 1.0)
    assert np.all(~np.isnan(ld_arr))

    # sanity checks
    assert ld_arr.shape[0] == ld_arr.shape[1]
    if ld_arr.shape[0] != ld_snps.shape[0]:
        raise ValueError('LD matrix has a different number of SNPs than the SNPs file')

    # get locus_allele
    ld_snps['locus_allele'] = ld_snps[['chromosome', 'position', 'allele1', 'allele2']].astype(str).agg('_'.join, axis=1)

    return ld_snps, ld_arr

#%%
def get_ld_for_snps(snps_to_slct, ld_snps, ld_arr):
    '''
        `snps_to_slct`: A list contains SNPs in the format of locus allele, e.g. chr_pos_ref_alt
        `ld_snps`: A pandas DataFrame containing locus_allele column in the same format as snps_to_slct
        `ld_arr`: A numpy ndarray

        Return:
            `locus_alleles_with_ld`: A list contains SNPs which have LD
            `ld_for_snps`: A LD matrix for SNPs which have LD

        Example:
            ld_prefix = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_15000001_18000001'
            ld_snps, ld_arr = load_alkes_ld(ld_prefix)

            snps_to_slct = ['18_15000071_C_G', '18_15000199_C_A', '18_15000609_G_T']

            snps_with_ld, ld_for_snps = get_ld_for_snps(snps_to_slct, ld_snps, ld_arr)
            print(ld_for_snps)
    '''
    matches = ld_snps['locus_allele'].isin(snps_to_slct).to_numpy()
    idxs = np.argwhere(matches).flatten().tolist()
    ld_for_snps = ld_arr[idxs, :][:, idxs]

    # check if all SNPs in GWAS summary has LD, if not, remove SNPs which don't have LD
    snps_with_ld = ld_snps.iloc[idxs]
    num_snps_with_ld = snps_with_ld.shape[0]

    if num_snps_with_ld != len(snps_to_slct):
        num_snps_removed = len(snps_to_slct) - num_snps_with_ld

        logging.info(f'Warning. Removed {num_snps_removed} SNPs since they do not have LD')

    locus_alleles_with_ld = snps_with_ld['locus_allele'].tolist()

    return locus_alleles_with_ld, ld_for_snps

#%%
def concat_ld_matrix(ld_snps_id1: list, ld_arr1: np.ndarray, ld_snps_id2: list, ld_arr2: np.ndarray):
    '''
        ld_snps_id1, ld_snps_id2: A list of SNP ids, must in physical order in genome
        ld_arr1, ld_arr2: A numpy ndarray, both with shape [n,n], [m,m]
        
        Example:
            ld_prefix1 = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_15000001_18000001'
            ld_prefix2 = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_18000001_21000001'

            ld_snps1, ld_arr1 = load_alkes_ld(ld_prefix1)
            ld_snps2, ld_arr2 = load_alkes_ld(ld_prefix2)

            ld_snps_id1 = ld_snps1['locus_allele'].tolist()
            ld_snps_id2 = ld_snps2['locus_allele'].tolist()

            final_ld_snps, final_ld_arr = concat_ld_matrix(ld_snps_id1, ld_arr1, ld_snps_id2, ld_arr2)
    '''
    t0 = time.time()

    last_snp_in_first_ld_region = ld_snps_id1[-1]
    first_snp_in_second_ld_region = ld_snps_id2[0]

    if last_snp_in_first_ld_region in ld_snps_id2 and first_snp_in_second_ld_region in ld_snps_id1:
        print('Found overlapping SNPs in two LD region')
        # find index of intersection SNP of both arrays
        arr2_idx_in_arr1 = ld_snps_id1.index(ld_snps_id2[0])
        arr1_idx_in_arr2 = ld_snps_id2.index(ld_snps_id1[-1])

        # generate lower left full zero matrix and upper right full zero matrix
        arr1_left_zero_num = arr2_idx_in_arr1
        arr2_left_zero_num = ld_arr2.shape[0]-arr1_idx_in_arr2-1

        # 这两个array互为转置
        lower_left_zero_arr = np.zeros((arr2_left_zero_num, arr1_left_zero_num))
        upper_right_zero_arr = np.zeros((arr1_left_zero_num, arr2_left_zero_num))

        # get left array and right array
        lower_left_arr = np.concatenate((lower_left_zero_arr, ld_arr2[(arr1_idx_in_arr2+1):ld_arr2.shape[0], 0:(arr1_idx_in_arr2+1)]), axis=1)
        right_arr = np.concatenate((upper_right_zero_arr, ld_arr2[:, (arr1_idx_in_arr2+1):ld_arr2.shape[0]]), axis=0)
        left_arr = np.concatenate((ld_arr1, lower_left_arr), axis=0)

        # get final array and SNP ids
        final_arr = np.concatenate((left_arr, right_arr), axis=1)
        final_ids = ld_snps_id1 + ld_snps_id2[(arr1_idx_in_arr2+1):len(ld_snps_id2)]

    else:
        print('No overlapping SNPs in two LD region. Assuming SNPs in two LD region are independent.')
        
        num_snps_in_region1 = ld_arr1.shape[0]
        num_snps_in_region2 = ld_arr2.shape[0]

        # generate lower left full zero matrix and upper right full zero matrix
        lower_left_zero_arr = np.zeros((num_snps_in_region2, num_snps_in_region1))
        upper_right_zero_arr = np.zeros((num_snps_in_region1, num_snps_in_region2))

        # concat to left ld and right ld
        left_arr = np.concatenate((ld_arr1, lower_left_zero_arr), axis=0)
        right_arr = np.concatenate((upper_right_zero_arr, ld_arr2), axis=0)

        # concat left ld and right ld to final ld
        final_arr = np.concatenate((left_arr, right_arr), axis=1)
        final_ids = ld_snps_id1 + ld_snps_id2

    print('Concatenating is done in %0.2f seconds'%(time.time() - t0))
    return final_ids, final_arr

#%%
def concat_ld_matrices(ld_prefixes: list):
    '''
        Main: Concatenating all LD matrices and SNPs defined in `ld_prefixes`.

        Arguments:
            ld_prefixes: A list contains LD prefixes which are used to be concatenated together. 
        
        Return: 
            concated_ld_snps: A pandas DataFrame object containing SNP IDs in the locus_allele column.
            concated_ld_matrix: A numpy ndarray object containing the concatenated LD matrix. 
    
        Example:
            ld_prefixes = [
                '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_15000001_18000001', 
                '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_18000001_21000001'
            ]
            final_ld_snps, final_ld_arr = concat_ld_matrices(ld_prefixes)
    '''
    num_total_ld_files = len(ld_prefixes)

    ## load first LD matrix and SNPs
    first_ld_prefix = ld_prefixes[0]
    ld_snps, ld_arr = load_alkes_ld(first_ld_prefix)

    concated_ld_matrix = ld_arr
    concated_ld_snp_ids = ld_snps['locus_allele'].tolist()

    ## concatenating the rest LD matrix to first LD matrix
    counter = 2
    for ld_prefix in ld_prefixes[1:]:
        logging.info(f'Concatenating {counter} / {num_total_ld_files} LD file')
        ld_snps, ld_arr_to_concat = load_alkes_ld(ld_prefix)
        snps_to_concat = ld_snps['locus_allele'].tolist()

        concated_ld_snp_ids, concated_ld_matrix = concat_ld_matrix(
            concated_ld_snp_ids, concated_ld_matrix, snps_to_concat, ld_arr_to_concat
        )

        counter += 1

    concated_ld_snps = pd.DataFrame({'locus_allele': concated_ld_snp_ids})

    print('All LD matrices are concatenated!')

    return concated_ld_snps, concated_ld_matrix

#%%
def example_codes():
    ### concat two LD regions using concat_ld_matrix
    ld_prefix1 = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_15000001_18000001'
    ld_prefix2 = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_18000001_21000001'

    ld_snps1, ld_arr1 = load_alkes_ld(ld_prefix1)
    ld_snps2, ld_arr2 = load_alkes_ld(ld_prefix2)

    ld_snps_id1 = ld_snps1['locus_allele'].tolist()
    ld_snps_id2 = ld_snps2['locus_allele'].tolist()

    final_ld_snps, final_ld_arr = concat_ld_matrix(ld_snps_id1, ld_arr1, ld_snps_id2, ld_arr2)

    ### concat multiple LD regions using concat_ld_matrices
    ld_prefixes = [
        '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_15000001_18000001', 
        '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_16000001_19000001', 
        '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_17000001_20000001', 
        '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/chr18_18000001_21000001'
    ]

    final_ld_snps, final_ld_arr = concat_ld_matrices(ld_prefixes)



#%%
def myfunc(args):
    input_file = args.input_file
    out_file = args.out_file

    pass

def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("--input_file", help="Project dir", required=True)

    parser.add_argument("--out_file", help="Project dir", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

def main():
    args = get_args()
    myfunc(args)

if __name__ == "__main__":
    main()
