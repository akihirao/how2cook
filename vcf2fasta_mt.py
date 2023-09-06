# vcf2fas_mt.py
# by HIRAO Akira
# A python script for converting vcf to fasta for mitogenome

'''
Usage:
    python vcf2fasta_mt -i sample.vcf (or sample.vcf.gz) -o output.fasta
    This script was tested under python 3.9 on conda with pyvcf
'''


import vcf
import argparse
import pandas as pd
import numpy as np


def format_str(x):
    '''
    Change numbers to a string, and change list to a string.
    '''
    if x is None:
        y = ''
    else:
        if isinstance(x, list):
            if not x:
                y = ';'.join(x)
            else:
                y = ''
        else:
            y = str(x)
    return y



def convert_gt(gt, ref_base, alt_base):
    '''
    conversion of genotype data format as Lumped  
        "0" -> ref_allele
        "1" -> alt_allele
        "." -> N
    '''
    if gt == ".":
        gt_Lumped =  "N"
    elif gt == '0':
        gt_Lumped = ref_base
    elif gt == '1':
        gt_Lumped = alt_base
    return gt_Lumped



def vcf2fasta_mt(inp,out):
    with open(out, 'w') as outfh:
        vcf_fh = vcf.Reader(filename=inp)
        samples_list = vcf_fh.samples
        no_sample = len(samples_list)
        loci_list = []
        tab_list = []
        df = pd.DataFrame([], index= samples_list)
        for snp_record in vcf_fh:
            ref_base = str(snp_record.REF)
            alt_base = str(snp_record.ALT)
            alt_base = alt_base[1]
            gt_list = []
            for sample in snp_record.samples:
                gt_raw = format_str(sample['GT'])
                gt_Lumped = convert_gt(gt_raw,ref_base, alt_base)
                gt_list.append(gt_Lumped)
            locus = format_str(snp_record.CHROM) + '_' + format_str(snp_record.POS)
            loci_list.append(locus)

            gt_list_df = pd.DataFrame(gt_list, index= samples_list)
            df = pd.concat([df, gt_list_df], axis = 1)

        no_loc = len(loci_list)
        space_list = np.repeat('  ',no_sample)
        sample_df = pd.DataFrame(samples_list,index=samples_list)
        space_df = pd.DataFrame(space_list,index=samples_list)
        df_out = pd.concat([sample_df, space_df, df], axis = 1)

        for index,item in df.iterrows():
            #print('>' + index)
            outfh.write('>' + index + '\n')
            sequence_list =item.to_list()
            sequence = ''.join(sequence_list)
            outfh.write(sequence + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Get Sample information.')
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    vcf2fasta_mt(args.input, args.output)
