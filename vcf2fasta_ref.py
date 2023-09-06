# vcf2fasta_ref.py
# by HIRAO Akira
# A python script for converting vcf to reference fasta

'''
Usage:
    python vcf2ref_fa.py -i sample.vcf (or sample.vcf.gz) -i output.fasta
    This script was tested under python 3.9 on conda with pyvcf
'''

import vcf
import argparse

def vcf2fasta_ref(inp, out):
    with open(out, 'w') as outfh:
        vcf_fh = vcf.Reader(filename=inp)
        samples_list = vcf_fh.samples
        no_sample = len(samples_list)
        ref_list = []
        for snp_record in vcf_fh:
            ref_base = str(snp_record.REF)
            ref_list.append(ref_base)

        sequence = ''.join(ref_list)
        outfh.write(">ref" + '\n')
        outfh.write(sequence + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Get Sample information.')
    parser.add_argument('-i', '--input', required = True)
    parser.add_argument('-o', '--output', required = True)
    args = parser.parse_args()
    vcf2fasta_ref(args.input, args.output)

