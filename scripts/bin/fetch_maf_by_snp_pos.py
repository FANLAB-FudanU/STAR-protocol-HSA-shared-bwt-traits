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

def get_alignments_from_maf(maf_file, ref_assembly):
    '''
        Position in MAF file is 0-based
    '''
    maf_file_reader = open_file(maf_file)
    for line in maf_file_reader:
        if not line.startswith('#'):
            if line.startswith('a'):
                alignment = {}
                alignment['score_line'] = line.strip()
                alignment['query_data'] = []
                alignment['tgt_data'] = []
            
            elif line.startswith(f's {ref_assembly}'):
                fields = line.strip().split()
                seq_name, start, aligned_size, strand, seq_len, aligned_seq = fields[1:]
                chrom = seq_name.split('.')[1]
                start = int(start) # 0-based
                aligned_size = int(aligned_size)
                seq_len = int(seq_len)
                end = start + aligned_size

                alignment['query_data'].extend(
                    [chrom, start, end, seq_name, seq_len, strand, aligned_seq]
                )
                
            elif line == '\n':
                yield alignment

            else:
                alignment['tgt_data'].append(line.strip().split())

def add_space_padding(
        seq_name, start, aligned_size, seq_len, 
        seq_name_padding_length = 28, start_padding_length = 10, 
        aligned_size_padding_length = 5, seq_len_padding_length = 10
    ):
    '''
    Example:
        seq_name = 'hg38.chr1'
        start = 12345569
        aligned_size = 10
        seq_len = 122229933

        padded_seq_name, padded_start, padded_aligned_seq, padded_seq_len = add_space_padding(seq_name, start, aligned_size, seq_len)
        print(padded_seq_name)
        print(padded_start)
        print(padded_aligned_seq)
        print(padded_seq_len)
    '''
    seq_name = str(seq_name).ljust(seq_name_padding_length)
    start = str(start).rjust(start_padding_length)
    aligned_size = str(aligned_size).rjust(aligned_size_padding_length)
    seq_len = str(seq_len).rjust(seq_len_padding_length)

    return seq_name, start, aligned_size, seq_len

def write_maf_alignments(
        alignments, out_file, 
        seq_name_padding_length = 28, start_padding_length = 10, 
        aligned_size_padding_length = 5, seq_len_padding_length = 10
    ):
    '''
    Args:
        alignments: MAF alignments loaded from MAF file by function `get_alignments_from_maf`
        out_file: Output file with `.maf` suffix.
    '''
    with open(out_file, 'w') as ofl:
        ofl.write('##maf version=1\n')

        for ali in alignments:
            score_line = ali['score_line']
            query_data = ali['query_data']
            tgt_data = ali['tgt_data']
            
            # write score line
            ofl.write(f'{score_line}\n')

            # write query data
            align_type = 's'
            chrom, start, end, seq_name, seq_len, strand, aligned_seq = query_data
            aligned_size = end - start
            padded_seq_name, padded_start, padded_aligned_seq, padded_seq_len = add_space_padding(
                seq_name, start, aligned_size, seq_len, 
                seq_name_padding_length, start_padding_length, aligned_size_padding_length, seq_len_padding_length
            )

            ofl.write(f'{align_type} {padded_seq_name} {padded_start} {padded_aligned_seq} {strand} {padded_seq_len} {aligned_seq}\n')

            # write tgt data
            for d in tgt_data:
                align_type = d[0]
                if align_type == 's' or align_type == 'e':
                    align_type, seq_name, start, aligned_size, strand, seq_len, aligned_seq = d
                    padded_seq_name, padded_start, padded_aligned_seq, padded_seq_len = add_space_padding(
                        seq_name, start, aligned_size, seq_len)

                    ofl.write(f'{align_type} {padded_seq_name} {padded_start} {padded_aligned_seq} {strand} {padded_seq_len} {aligned_seq}\n')
                else:
                    # 以非s开头的比对结果取出来是未比对上的序列，原样输出即可
                    seq_name = d[1]
                    padded_seq_name = seq_name.ljust(seq_name_padding_length)
                    ofl.write(f'{align_type} {padded_seq_name} ')
                    ofl.write('%s\n'%(' '.join(d[2:])))

            ofl.write('\n')

        ofl.write('#eof\n')


#%%
def extract_seq_from_alignment(tgt_start, alignment):
    res_maf = []
    res_maf.append(alignment['score_line'])

    chrom, start, end, seq_name, seq_len, strand, aligned_seq = alignment['query_data']

    # get genomic position of aligned bases
    aligned_seq_idxs = [] 
    curr_idx = start
    for s in aligned_seq:
        if s != '-':
            aligned_seq_idxs.append(curr_idx)
            curr_idx += 1
        else:
            aligned_seq_idxs.append('-')

    base_idx = aligned_seq_idxs.index(tgt_start)
    query_base = aligned_seq[base_idx]

    # add query line
    res_maf.append(['s', seq_name, tgt_start, 1, strand, seq_len, query_base])

    # add tgt lines of other species
    idx = 0
    last_line_is_unaligned = False
    for tgt_line in alignment['tgt_data']:
        if tgt_line[0] == 's':
            align_type, seq_name, start, aligned_size, strand, seq_len, aligned_seq = tgt_line
            start = int(start)
            
            tgt_base = aligned_seq[base_idx]

            # get start in tgt species
            aligned_seq_idxs = [] 
            curr_idx = start
            for s in aligned_seq:
                if s != '-':
                    aligned_seq_idxs.append(curr_idx)
                    curr_idx += 1
                else:
                    aligned_seq_idxs.append('-')
            
            new_start = aligned_seq_idxs[base_idx]

            if tgt_base == '-': # 当前位置无比对上的碱基, 因此align type为e，且下一行i开头的行不需要输出
                status = 'C'
                last_line_is_unaligned = True
                last_line_seq_name = seq_name
                res_maf.append(['e', seq_name, new_start, 0, strand, seq_len, status])

            else:
                res_maf.append(['s', seq_name, new_start, 1, strand, seq_len, tgt_base])
        
        else:
            # 上一行为s开头的比对结果取出来是未比对上的序列，则后面一行以i开头的行不输出
            if last_line_is_unaligned and tgt_line[0] == 'i':
                seq_name = tgt_line[1]
                if last_line_seq_name == seq_name:
                    last_line_is_unaligned = False
                    continue

            res_maf.append(tgt_line)

        idx += 1

    return res_maf

#%%
def fetch_maf_by_snp_pos(args):
    raw_maf_file = args.raw_maf_file
    snp_pos_file = args.snp_pos_file
    human_assembly = args.human_assembly
    out_prefix = args.out_prefix

    # read snp pos
    snp_poss = []
    print(f'Reading snp pos from {snp_pos_file}')
    with open(snp_pos_file, 'r') as ifl:
        for l in ifl:
            line = l.strip().split()
            chrom, start, end = line[0:3]
            start, end = int(start), int(end)
            snp_poss.append([chrom, start, end])

    num_total_snps = len(snp_poss)
    print(f'Total {num_total_snps} positions are loaded')

    alignments = get_alignments_from_maf(raw_maf_file, human_assembly)

    num_not_found_snps = 0
    curr_pos_idx = 0
    snp_poss_updated = [] # 新增一列，表示这个位置的序列是否取出来了
    fetched_maf = []
    found = False

    p_iter = 1000

    print(f'Extracting maf from {raw_maf_file}')
    for alignment in alignments:
        query_data = alignment['query_data']
        chrom, start, end, seq_name, seq_len, strand, aligned_seq = query_data

        # print(f'alignment_idx is {alignment_idx}')
        while curr_pos_idx < num_total_snps:
            # print(f'curr_pos_idx is {curr_pos_idx}')

            # get curr snp
            tgt_chrom, tgt_start, tgt_end = snp_poss[curr_pos_idx]

            if tgt_end <= end:
                # print('tgt_end and end is', tgt_end, end)

                # alignment contains this snp
                if (tgt_chrom == chrom) and (tgt_start >= start):
                    # get maf
                    if tgt_start == 1081816:
                        res = alignment
                    # print(f'tgt_start is {tgt_start}')
                    maf = extract_seq_from_alignment(tgt_start, alignment)

                    found = True
                    snp_poss_updated.append([tgt_chrom, tgt_start, tgt_end, found])
                    curr_pos_idx += 1

                    fetched_maf.append(maf)

                    if curr_pos_idx % p_iter == 0:
                        print(f'{curr_pos_idx} / {num_total_snps} SNPs are extracted')

                # alignment does not contain this snp, need to process next snp
                else:
                    found = False            
                    snp_poss_updated.append([tgt_chrom, tgt_start, tgt_end, found])
                    curr_pos_idx += 1
                    num_not_found_snps += 1

            # curr alignment is at upstream of this snp, process next alignment
            else:
                break
        
        if curr_pos_idx == num_total_snps:
            print('MAF of all snps are successfully extracted!')
            break
    
    fetched_maf_file = f'{out_prefix}.maf.tsv.gz'
    snp_pos_stature_file = f'{out_prefix}.snp_status.tsv'

    # write into file
    with gzip.open(fetched_maf_file, 'wt') as ofl:
        ofl.write('##maf version=1\n')
        for maf in fetched_maf:
            score_line = maf[0]
            ofl.write(f'{score_line}\n')
            for l in maf[1:]:
                ofl.write('%s\n'%('\t'.join([str(x) for x in l])))
            ofl.write('\n')

    # write snp_poss_updated
    with open(snp_pos_stature_file, 'w') as ofl:
        ofl.write('chrom\tstart\tend\tfound_maf\n')
        for l in snp_poss_updated:
            ofl.write('%s\n'%('\t'.join([str(x) for x in l])))

    print(f'Done. maf are written into {fetched_maf_file}')
    print(f'Done. status of position are written into {snp_pos_stature_file}')

#%%
def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("--raw_maf_file", help="MAF file", required=True)
    parser.add_argument("--snp_pos_file", help="Bed file containing positions of SNPs. Only first three columns are used. Additional columns are allowed but not used", required=True)
    parser.add_argument("--human_assembly", help="Project dir", required=True)
    parser.add_argument("--out_prefix", help="Project dir", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

def main():
    args = get_args()
    fetch_maf_by_snp_pos(args)

if __name__ == "__main__":
    main()
