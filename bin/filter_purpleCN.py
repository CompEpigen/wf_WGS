#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='Copy number file')
parser.add_argument('-o', type = str, help='Filtered copy number file')
parser.add_argument('--vcf', type = str, help='VCF of structural variants. Short CNAs are kept only if they are supported by the filtered SVs.')
args = parser.parse_args()

df = pd.read_csv(args.i,sep='\t')

breakends = {}
reader = vcfpy.Reader.from_path(args.vcf)
for record in reader:
    chr= str(record.CHROM)
    pos = record.POS
    chr2 = str(record.ALT[0].mate_chrom)
    pos2 = record.ALT[0].mate_pos
    if not chr in breakends: breakends[chr] = set()
    if not chr2 in breakends: breakends[chr2] = set()
    breakends[chr].add(pos)
    breakends[chr2].add(pos2)

def pos_in_breakends(chr,pos,breakends):
    if str(chr) in breakends:
        for x in breakends[str(chr)]:
            if abs(x-pos)<8: return True
    return False


i=0
while i <df.shape[0]-1:
    if df.loc[i,"chromosome"]==df.loc[i+1,"chromosome"]:
        if ( abs(df.loc[i,"copyNumber"]-df.loc[i+1,"copyNumber"])<0.15 or (abs(df.loc[i,"copyNumber"]-df.loc[i+1,"copyNumber"])<0.25 and df.loc[i+1,"bafCount"]<40) ) \
        and ( abs(df.loc[i,"baf"]-df.loc[i+1,"baf"])<0.06 or df.loc[i+1,"bafCount"]<50 ): # same chromosome, same copy number, no LOH
            f1 = df.loc[i,"end"] -  df.loc[i,"start"]
            f2 = df.loc[i+1,"end"] -  df.loc[i+1,"start"]
            df.loc[i,"copyNumber"] = (df.loc[i,"copyNumber"] * f1 + df.loc[i+1,"copyNumber"] * f2) / (f1+f2) # average copy numbers of the two merged segments
            df.loc[i,"end"] = df.loc[i+1,"end"]
            df.drop(i+1,axis=0,inplace=True)
            df.reset_index(drop=True,inplace=True)
        elif df.loc[i+1,"method"]=="UNKNOWN" or df.loc[i+1,"end"]-df.loc[i+1,"start"]<400:
            df.drop(i+1,axis=0,inplace=True)
            df.reset_index(drop=True,inplace=True)
        elif df.loc[i+1,"end"]-df.loc[i+1,"start"]<20000 and \
            ((not pos_in_breakends(df.loc[i+1,"chromosome"],df.loc[i+1,"start"],breakends)) or (not pos_in_breakends(df.loc[i+1,"chromosome"],df.loc[i+1,"end"],breakends))):
            df.loc[i,"end"] = df.loc[i+1,"end"]
            df.drop(i+1,axis=0,inplace=True)
            df.reset_index(drop=True,inplace=True)
        elif (not pos_in_breakends(df.loc[i+1,"chromosome"],df.loc[i+1,"start"],breakends)) and abs(df.loc[i,"baf"]-df.loc[i+1,"baf"])<0.1 and abs(df.loc[i,"copyNumber"]-df.loc[i+1,"copyNumber"])<0.07:
            df.loc[i,"end"] = df.loc[i+1,"end"]
            df.drop(i+1,axis=0,inplace=True)
            df.reset_index(drop=True,inplace=True)
        else:
            i+=1
    else:
        i+=1

df.to_csv(args.o,sep="\t",index=False)