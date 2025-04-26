#!/usr/bin/env python
#%%
import os
import sys
import gzip
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """
        可以将参数的默认选项打印在帮助文档中
    """
    pass

def check_files(*args):
    for afile in args:
        if not os.path.exists(afile):
            raise Exception("Error! File `%s` not exists! Please check it!"%(afile))

def get_blocks_from_maf(maf_file, human_assembly, assembly_species_table):
    with gzip.open(maf_file, 'rt') as ifl:
        for line in ifl:
            if not line.startswith('#'):
                if line.startswith('a'):
                    curr_block = []

                if line == '\n':
                    yield curr_block
                else:
                    if line.startswith('s') or line.startswith('e'): # only keep s and e line
                        line = line.strip().split()
                        
                        # only keep species that in assembly_species_table
                        assembly_version = line[1].split('.')[0]
                        if (assembly_version == human_assembly or 
                            assembly_version in assembly_species_table):
                            curr_block.append(line)


def reformat_maf(args):
    input_file = args.input_file
    species_assembly_file = args.species_assembly_file
    out_file = args.out_file
    human_assembly = args.human_assembly

    print(f'input_file is {input_file}')
    print(f'species_assembly_file is {species_assembly_file}')
    print(f'out_file is {out_file}')

    # read species_assembly_file
    assembly_species_table = {}
    with open(species_assembly_file, 'r') as ifl:
        ifl.readline()
        for l in ifl.readlines():
            species, assembly_version = l.strip().split('\t')
            assembly_species_table[assembly_version] = species

    print(f'Total {len(assembly_species_table)} species are loaded')

    # reformat results
    print('INFO | Reformatting blocks')
    blocks_iter = get_blocks_from_maf(input_file, human_assembly, assembly_species_table)

    out_file_writer = gzip.open(out_file, 'wt')
    out_file_writer.write('human.chrom\thuman.pos\thuman.seq')
    for species_name in assembly_species_table.values():
        out_file_writer.write(f'\t{species_name}.linetype\t{species_name}.chrom\t{species_name}.pos\t{species_name}.seq')
    out_file_writer.write('\n')

    # header: human.chrom, human.pos, human.seq, 
    # {species}.linetype, {species}.chrom, {species}.pos, {species}.seq
    p_iter = 10000
    counter = 0
    for block in blocks_iter:
        refmtted_block = {}
        for line in block:
            assembly_version = line[1].split('.')[0]
            if assembly_version == human_assembly:
                human_chrom = line[1].split('.')[1]
                human_pos = line[2] # 0-seqd
                human_seq = line[-1].upper()[0]
                refmtted_block['human'] = [human_chrom, human_pos, human_seq]
            else:
                species_name = assembly_species_table[assembly_version]
                species_linetype = line[0]
                species_chrom = line[1].split('.')[1]
                species_pos = line[2]
                species_seq = line[-1].upper()[0]
                refmtted_block[species_name] = [species_linetype, species_chrom, species_pos, species_seq]
        
        # write into file
        chrom, pos, seq = refmtted_block['human']
        out_file_writer.write(f'{chrom}\t{pos}\t{seq}')

        # write info of other species
        for species_name in assembly_species_table.values():
            if species_name in refmtted_block:
                linetype, chrom, pos, seq = refmtted_block[species_name]
            else:
                linetype, chrom, pos, seq = 'NA', 'NA', 'NA', 'NA'
            out_file_writer.write(f'\t{linetype}\t{chrom}\t{pos}\t{seq}')
        out_file_writer.write('\n')
        out_file_writer.flush()

        counter += 1
        if counter % p_iter == 0:
            print(f"{counter} blocks are reformatted")

    print(f'INFO | Done. out file is {out_file}')


#%%
def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("-i", "--input_file", help="Gzipped maf file obtained from UCSC table browser", required=True)
    parser.add_argument("--species_assembly_file", help="Header is required", required=True)
    parser.add_argument("--human_assembly", help="Assembly version of human genome. Can be hg19 or hg38", required=True, default='hg19')
    parser.add_argument("--out_file", help="Output file will be gzipped. Position is 0-based", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    reformat_maf(args)

    return args

def main():
    args = get_args()

if __name__ == "__main__":
    main()
