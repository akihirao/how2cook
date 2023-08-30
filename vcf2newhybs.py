# vcf2newhybs.py
# by HIRAO Akira
# A python script for converting vcf to an input file for NewHybrids (Anderson and Thompson 2002)

'''
Usage:
    python vcf2newhybs.py -i sample.vcf (or sample.vcf.gz) -o output.txt
    This script was tested under python 3.9 on conda with pyvcf
'''


import vcf
import random
import argparse
import pandas as pd


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



def convert_gt(gt):
    '''
    conversion of genotype data format as Lumped
        "0/0" or "0|0" -> "11"
        "0/1" or "0|1" -> "12"
        "1/1" or "1|1" -> "22"
        "./.",".|." or "." -> 0
    '''
    if gt[0] == ".":
        gt_Lumped = 0
    elif gt == '0/0' or gt == '0|0':
        gt_Lumped = 11
    elif gt == '0/1' or gt == '0|1':
        gt_Lumped = 12
    elif gt == '1/1' or gt == '1|1':
        gt_Lumped = 22
    return gt_Lumped


def prior_def_species(sample):
    '''
    prior definition of pure species as newhybrids format: z0 or z1
    This function is specific for labeling snow crab or red snow crab!!
    '''
    if sample[0:2] == "CJ":
        z_value = "z1"
    else:
        z_value = "z0"
    return z_value



def vcf2newhybs(inp, out):
    no_selected_loci = 300
    with open(out, 'w') as outfh:
        vcf_fh = vcf.Reader(filename=inp)
        samples_list = vcf_fh.samples
        z_list = []
        for each_sample in samples_list:
            z_value = prior_def_species(each_sample)
            z_list.append(z_value)
        no_sample = len(samples_list)
        loci_list = []
        df = pd.DataFrame([], index= samples_list)
        for snp_record in vcf_fh:
            gt_list = []
            for sample in snp_record.samples:
                gt_Lumped = convert_gt(sample['GT'])
                gt_list.append(gt_Lumped)
            locus = format_str(snp_record.CHROM) + '_' + format_str(snp_record.POS)
            loci_list.append(locus)
            gt_list_df = pd.DataFrame(gt_list, index= samples_list)
            df = pd.concat([df, gt_list_df], axis = 1)

        # randomly selection of xxx loci
        no_loci = len(loci_list)
        selected_loci_id = sorted(random.sample(list(range(0,no_loci)),no_selected_loci))
        df_selected = df.iloc[:,selected_loci_id]
        loci_selected_list = []
        for extract_id in selected_loci_id:
            selected_locus = loci_list[extract_id]
            loci_selected_list.append(selected_locus)

        outfh.write('NumIndivs' + '\t' + str(no_sample) + '\n')
        outfh.write('NumLoci' + '\t' + str(no_selected_loci) + '\n')
        outfh.write('Digits' + '\t' + '1' + '\n')
        outfh.write('Format Lumped' + '\n')
        outfh.write('\n')
        outfh.write('LocusNames'+'\t')
        outfh.write('\t'.join(loci_selected_list) + '\n')
        outfh.write('\n')

    id_df = pd.DataFrame(list(range(1, no_sample + 1)),index=samples_list)
    n_df = pd.DataFrame(['n' for i in range(no_sample)],index=samples_list)
    sample_df = pd.DataFrame(samples_list,index=samples_list)
    z_df = pd.DataFrame(z_list, index=samples_list)
    df_out = pd.concat([id_df, n_df, sample_df, z_df, df_selected], axis = 1)
    df_out.to_csv(out, mode='a', sep = '\t', index=False, header=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Get Sample information.')
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    vcf2newhybs(args.input, args.output)

