#!/usr/bin/env python
#%%
import os
import sys
import time
import gzip
import argparse
from collections import OrderedDict

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """
        可以将参数的默认选项打印在帮助文档中
    """
    pass

def check_files(*args):
    for afile in args:
        if not os.path.exists(afile):
            raise Exception("Error! File `%s` not exists! Please check it!"%(afile))

def timer(func):
    '''
        Run function and report used time
    '''
    def wrapper():
        start_time = time.time()

        func()

        end_time = time.time()
        used_time = end_time - start_time
        used_mins = used_time / 60
        print('Total time: %.8s seconds'%(used_time))
    return wrapper

def open_file(afile):
    if afile.endswith('gz'):
        file_reader = gzip.open(afile, 'rt')
    else:
        file_reader = open(afile)
    return file_reader

def find_human_specific_allele(human_ref, human_alt, other_species_linetype, other_species_seq):
    '''
        other_species_linetype: Referred to linetype of MAF file
        other_species_seq: Referred to aligned sequence of MAF file

        Example: 
            human_ref = 'A'
            human_alt = 'G'
            other_species_linetype = 's'
            other_species_seq = 'G'
            find_human_specific_allele(human_ref, human_alt, other_species_linetype, other_species_seq)
    '''

    human_alleles = [human_ref, human_alt]

    # other_species_linetype is e or i, or NA
    if ( other_species_linetype == 'e' or 
        other_species_linetype == 'i' or 
        other_species_linetype == 'NA' ):
        other_species_alleles = []

    elif other_species_linetype == 's':
        other_species_alleles = [other_species_seq]

    else:
        raise Exception(f'Error! Not valid linetype {other_species_linetype}. Can only be one of s, e, i or NA')

    # remove alleles of other species
    human_specific_alleles = list(sorted(set(human_alleles) - set(other_species_alleles)))

    return human_specific_alleles


#%%
def find_hsa_snps(args):
    '''
        Find SNPs with human specific allele (hsa)
    '''
    input_file = args.input_file
    species_assembly_file = args.species_assembly_file
    out_file = args.out_file

    print(f'input_file is {input_file}')
    print(f'species_assembly_file is {species_assembly_file}')
    print(f'out_file is {out_file}')

    # read species_assembly_file
    species_list = []
    with open(species_assembly_file, 'r') as ifl:
        ifl.readline()
        for l in ifl.readlines():
            species_name, assembly_version = l.strip().split('\t')
            species_list.append(species_name)

    print(f'Total {len(species_list)} species are loaded')

    out_file_writer = gzip.open(out_file, 'wt')
    out_file_writer.write('chrom\tstart\tend\tref\talt\t')
    # out_file_writer.write('chrom\tstart\tend\tref\talt\tchimp.seq\tgorilla.seq\torangutan.seq\thuman_specific_alleles\n')
    for species_name in species_list:
        out_file_writer.write(f'{species_name}.linetype\t{species_name}.seq\t')
    out_file_writer.write('human_specific_alleles\n')

    counter = 0
    p_iter = 10000
    species_idxs = OrderedDict()

    file_reader = open_file(input_file)
    header = file_reader.readline().strip().split('\t')

    # find index of each species
    for species_name in species_list:
        start_idx = header.index(f'{species_name}.linetype')
        species_idxs[species_name] = start_idx

    # find hsa for each snp
    for l in file_reader:
        line = l.strip().split()
        chrom, start, end, ref, alt, human_seq = line[0:6]

        # find hsa
        if human_seq != 'NA': 
            hsas = set(['A', 'C', 'G', 'T'])
            for species_name, start_idx in species_idxs.items():
                species_linetype, species_seq = line[start_idx], line[start_idx + 3]
                species_hsa = find_human_specific_allele(ref, alt, species_linetype, species_seq)
                hsas &= set(species_hsa)
            hsas = list(sorted(hsas))
        else: # 该snp在maf文件中不存在对应的比对结果
            hsas = ['NA']

        if not hsas:
            hsas = ['NA']

        out_file_writer.write(f'{chrom}\t{start}\t{end}\t{ref}\t{alt}\t')
        for species_name, start_idx in species_idxs.items():
            species_linetype, species_seq = line[start_idx], line[start_idx + 3]
            out_file_writer.write(f'{species_linetype}\t{species_seq}\t')

        # out_file_writer.write(f'{chimp_seq}\t{gorilla_seq}\t{orangutan_seq}\t')
        out_file_writer.write(f'{",".join(hsas)}\n')
        out_file_writer.flush()

        counter += 1
        if counter % p_iter == 0:
            print(f'{counter} SNPs are processed')

    file_reader.close()
    out_file_writer.close()

    print(f'Done. output file is {out_file}')

    # output columns: chrom, start, end, ref and alt
    # chimp.seq, gorilla.seq, orangutan.seq, human_specific_allele


def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("--input_file", help="Merged human and non-human-primates alleles. Ref and alt allele of human are required. First five columns are: chrom, start, end, ref and alt.", required=True)
    parser.add_argument("--species_assembly_file", help="Species assembly lookup file. Header is required", required=True)
    parser.add_argument("--out_file", help="Output file", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    find_hsa_snps(args)

    return args

@timer
def main():
    args = get_args()

if __name__ == "__main__":
    main()

#%%%
sys.exit()

#%%
