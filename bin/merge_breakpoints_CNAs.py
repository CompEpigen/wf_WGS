#!/usr/bin/env python3


import os
import sys
import numpy as np
import pandas as pd
import vcfpy

def get_breakpoint_info(record):
    chr = record.CHROM
    pos = record.POS
    if record.INFO["SVTYPE"] in ["DEL","DUP","INS"]:
        chr2 = chr
        pos2 = record.INFO["END"]
        if record.INFO["SVTYPE"]=="DEL":
            orientation = "-"
            orientation2 = "+"
        elif record.INFO["SVTYPE"]=="DUP":
            orientation = "+"
            orientation2 = "-"
        else:
            orientation = "-"
            orientation2 = "+"
    else:
        chr2 = record.ALT[0].mate_chrom
        pos2 = record.ALT[0].mate_pos
        orientation = record.ALT[0].orientation
        orientation2 = record.ALT[0].mate_orientation
    return ( (chr,pos,orientation) , (chr2,pos2,orientation2) )

d={"sample":[],"chr1":[],"pos1":[],"orientation1":[],"chr2":[],"pos2":[],"orientation2":[]}
l=[]
for x in sorted(os.listdir(sys.argv[1])):
    patient = x.split(".")[0]
    reader = vcfpy.Reader.from_path(os.path.join(sys.argv[1],x))
    for record in reader:
        d["sample"].append(patient)
        ( (chr,pos,orientation) , (chr2,pos2,orientation2) ) = get_breakpoint_info(record)
        d["chr1"].append(chr)
        d["pos1"].append(pos)
        d["orientation1"].append(orientation)
        d["chr2"].append(chr2)
        d["pos2"].append(pos2)
        d["orientation2"].append(orientation2)
df = pd.DataFrame(d)
df.to_csv("breakpoints.tsv",sep="\t",index=False)

d_CNA={"sample":[],"chr":[],"start":[],"end":[],"cn":[]}
for x in os.listdir(sys.argv[2]):
    patient = x.split(".")[0]
    df_CNA = pd.read_csv(os.path.join(sys.argv[2],x),sep="\t")
    for i in df_CNA.index:
        chr=df_CNA.loc[i,"chromosome"]
        start=df_CNA.loc[i,"start"]
        end=df_CNA.loc[i,"end"]
        cn=df_CNA.loc[i,"copyNumber"]
        if cn<2.5 and cn>1.5: continue
        cn = round(cn)
        if end-start>300:
            d_CNA["sample"].append(patient)
            d_CNA["chr"].append(chr)
            d_CNA["start"].append(start)
            d_CNA["end"].append(end)
            d_CNA["cn"].append(cn)



df_CNA=pd.DataFrame(d_CNA)
df_CNA.to_csv("CNAs.tsv",sep="\t",index=False)
