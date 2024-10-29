#!/usr/bin/env python3

import os
import argparse
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Input VCF file')
parser.add_argument('-o', type = str, help='Output vcf file')
args = parser.parse_args()

reader = vcfpy.Reader.from_path(args.i)

writer = vcfpy.Writer.from_path(args.o,reader.header)
for record in reader:
    chr = record.CHROM
    pos = record.POS
    #print(record)
    #print(record.calls[0])
    AD = record.calls[0].data["AD"]
    if len(AD)<2: continue
    if AD[0]>5 and AD[1] > 5 and ( AD[0] / (AD[0]+AD[1]) >= 0.30) and (AD[0] / (AD[0]+AD[1]) <= 0.70):
        if len(record.REF)==1 and len(record.ALT)==1 and len(record.ALT[0].value)==1: # Only keep SNVs
            writer.write_record(record)