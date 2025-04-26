#!/usr/bin/env python

import os
import sys
import gzip
import argparse
import subprocess
from operator import itemgetter

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

def convert_genomic_versions(args):
    ## inputs
    input_file = args.input_file
    coord_lookup_table_file = args.coord_lookup_table_file
    out_file = args.out_file

    ## cmds
    coord_lookup_table_reader = open_file(coord_lookup_table_file)
    coord_lookup_table_header = coord_lookup_table_reader.readline().strip().split()

    chrom_idx = coord_lookup_table_header.index('chrom')
    hg19_pos_idx = coord_lookup_table_header.index('hg19_pos')
    hg19_ref_idx = coord_lookup_table_header.index('hg19_ref')
    hg38_pos_idx = coord_lookup_table_header.index('hg38_pos')

    ## read coord_lookup_table
    coord_lookup_table = {}
    for l in coord_lookup_table_reader:
        line = l.strip().split()
        chrom = line[chrom_idx]
        hg19_pos, hg19_ref, hg38_pos = line[hg19_pos_idx], line[hg19_ref_idx], line[hg38_pos_idx]

        hg38_id = f'{chrom}_{hg38_pos}'
        coord_lookup_table[hg38_id] = [hg19_pos, hg19_ref]

    ## convert coords
    input_file_reader = open_file(input_file)
    input_file_header = input_file_reader.readline()

    counter = 0
    p_iter = 100000
    num_success = 0
    converted_snps = []

    for l in input_file_reader:
        line = l.strip().split()

        chrom, hg38_start, hg38_end, ref, alt, hsa = line
        hg38_id = f'{chrom}_{hg38_start}'
        if hg38_id in coord_lookup_table:
            hg19_pos, hg19_ref = coord_lookup_table[hg38_id]

            if hg19_ref == ref:
                hg19_alt = alt
            else:
                hg19_alt = ref

            converted_snps.append([chrom, int(hg19_pos), int(hg19_pos), hg19_ref, hg19_alt, hsa])
            num_success += 1

        else:
            print(f'Warning! SNP {hg38_id} is not in coord lookup table!')
            continue
        
        counter += 1
        if counter % p_iter == 0:
            print(f'{counter} SNPs done')

    print(f'Total {num_success} out of {counter} SNPs are converted successfully!')
    
    ## sort results
    print('Sorting results')
    converted_snps_sorted = list(sorted(converted_snps, key=itemgetter(0, 1, 2)))

    ## write results
    print('Writing results')
    with open(out_file, 'w') as ofl:
        ofl.write(input_file_header)
        for l in converted_snps_sorted:
            ofl.write('%s\n'%('\t'.join([str(x) for x in l])))

    ## compress and create index for out_file
    print('Compressing results and create index')
    cmd = f'bgzip -f {out_file}; tabix -S1 -s1 -b2 -e3 {out_file}.gz'
    print(f'cmd is {cmd}')
    run_shell_cmd(cmd)
    
    print(f'Done. out file is {out_file}.gz')

def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter, description='Convert hg38 position to hg19 position for SNPs with HSAs file')
    parser.add_argument("--input_file", help="Columns of input file are: chrom, start, end, ref, alt, human_specific_alleles. Should be output file of scripts/human_specific_alleles/get_snps_with_human_specific_alleles.sh", required=True)
    parser.add_argument("--coord_lookup_table_file", help="Coord lookup table used to convert hg38 pos to hg19 pos. Must contain at least three columns: chrom, hg19_pos and hg38_pos", required=True)
    parser.add_argument("--out_file", help="Output file. No gz suffix. It will be compressed and indexed automatically", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

def main():
    args = get_args()
    convert_genomic_versions(args)

if __name__ == "__main__":
    main()
