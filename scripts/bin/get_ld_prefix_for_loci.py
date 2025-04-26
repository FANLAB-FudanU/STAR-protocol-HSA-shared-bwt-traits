#!/usr/bin/env python
#%%
import os
import sys
import gzip
import argparse
import subprocess

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """
        可以将参数的默认选项打印在帮助文档中
    """
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

#%%
# test
# os.chdir('/public/home/fan_lab/wangjie/modern_human_specific_allele')

# ld_dir = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/'

# loci_file = 'tmp/loci.tsv'


#%%
def find_ld_prefix(ld_dir: str, chrom: str, min_pos: int):
    '''
        ld_dir = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/'
        chrom='12'
        min_pos = 40590585; 
        find_ld_prefix(ld_dir, chrom, min_pos)
    '''

    bp_starts = list(range(1, 252000001, 1000000))
    bp_ends = [x+3000000 for x in bp_starts]
    i = max([i for i,x in enumerate(bp_starts) if x<=min_pos])
    
    prefix = "chr%s_%d_%d" % (chrom, bp_starts[i], bp_ends[i])
    ld_prefix = os.path.join(ld_dir, prefix)

    return ld_prefix

#%%
def find_ld_prefix_v2(ld_dir: str, chrom: str, start: int, stop: int):
    '''
        ld_dir = '/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/'
        chrom='12'
        start = 91860962
        stop = 97339301
        find_ld_prefix_v2(ld_dir, chrom, start, stop)
    '''
    ld_file_size = 3000000 # 3 Mb

    bp_starts = list(range(1, 252000001, 1000000))
    bp_ends = [x+ld_file_size for x in bp_starts]
    i = max([i for i,x in enumerate(bp_starts) if x<=start])

    ld_start = bp_starts[i]
    ld_stop = bp_ends[i]

    prefixes = []
    prefix = "chr%s_%d_%d" % (chrom, ld_start, ld_stop)
    prefixes.append(prefix)

    # This locus spans at least two LD files
    while ld_stop < stop:
        ld_start = ld_start + 1000000
        ld_stop = ld_stop + 1000000

        prefix = "chr%s_%d_%d" % (chrom, ld_start, ld_stop)
        prefixes.append(prefix)

    ld_prefixes = [os.path.join(ld_dir, prefix) for prefix in prefixes]

    return ld_prefixes

#%%
def ld_files_exist(ld_prefixes):
    def _ld_exist(ld_prefix):
        snp_file = f'{ld_prefix}.gz'
        ld_file = f'{ld_prefix}.npz'
        ld_file2 = f'{ld_prefix}.npz2'

        if os.path.exists(snp_file):

            if os.path.exists(ld_file):
                file_existed = 'npz_file'
            elif os.path.exists(ld_file2):
                file_existed = 'npz2_file'
            else:
                file_existed = 'NA'

        else:
            file_existed = 'NA'
        
        return file_existed

    file_existed = list(set([_ld_exist(ld_prefix) for ld_prefix in ld_prefixes]))

    return ';'.join(file_existed)

#%%
def get_ld_prefix_for_loci(args):
    ld_dir = args.ld_dir
    loci_file = args.loci_file
    chrom_col_name = args.chrom_col_name
    start_col_name = args.start_col_name
    stop_col_name = args.stop_col_name
    out_file = args.out_file

    print(f'ld_dir is {ld_dir}')
    print(f'loci_file is {loci_file}')
    print(f'chrom_col_name is {chrom_col_name}')
    print(f'start_col_name is {start_col_name}')
    print(f'stop_col_name is {stop_col_name}')
    print(f'out_file is {out_file}')

    out_file_writer = open(out_file, 'w')

    ## read loci and find ld file
    print('Finding LD prefix for loci')
    counter = 0
    p_iter = 100
    with open(loci_file, 'r') as ifl:
        # get header
        header = ifl.readline().strip().split()

        # get chrom_idx
        if chrom_col_name in header:
            chrom_idx = header.index(chrom_col_name)
        else:
            if 'chr' in header:
                chrom_idx = header.index('chr')
            else:
                raise Exception(f'Error! Input loci file does not contain {chrom_col_name} column')
        
        # get start and stop idx
        start_idx = header.index(start_col_name)
        stop_idx = header.index(stop_col_name)

        # write header
        out_file_writer.write('%s\tnum_LD_regions\tld_files_status\tld_prefix\n'%('\t'.join(header)))

        # check ld files status
        for l in ifl:
            line = l.strip().split()
            chrom, start, stop = line[chrom_idx], line[start_idx], line[stop_idx]

            ld_prefixes = find_ld_prefix_v2(ld_dir, chrom, int(start), int(stop))
            num_ld_regions = len(ld_prefixes)

            ld_files_status = ld_files_exist(ld_prefixes)
            ld_prefixes_str = ';'.join(ld_prefixes)

            out_file_writer.write('%s\t%s\t%s\t%s\n'%(
                '\t'.join(line), num_ld_regions, ld_files_status, ld_prefixes_str
            ))
            out_file_writer.flush()

            counter += 1
            if counter % p_iter == 0:
                print(f'{counter} loci are done')

    print(f'Done. out file is {out_file}')


def get_args():
    ## get ld prefix from LD produced by AlkesGroup
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("--ld-dir", help="Directory containing LD files", required=True, default='/public/home/fan_lab/wangjie/PublicDB/GWAS/AlkesGroup_UKB/LD/01.downloaded/raw_data/')
    parser.add_argument("--loci-file", help="Locus file. chrom(chr) and start columns are required. No `chr` string in chrom column. If other colnames are used, please specify them with --chrom-col-name and --start-col-name parameters", required=True)
    parser.add_argument("--chrom-col-name", help="Colname of chromosome of each locus", required=False, default='chrom')
    parser.add_argument("--start-col-name", help="Colname of start position of each locus", required=False, default='start')
    parser.add_argument("--stop-col-name", help="Colname of stop position of each locus", required=False, default='stop')

    parser.add_argument("--out-file", help="Output file with LD status and ld prefix added. Locus with multiple LD regions are separated by semicolon", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

def main():
    args = get_args()
    get_ld_prefix_for_loci(args)

if __name__ == "__main__":
    main()
