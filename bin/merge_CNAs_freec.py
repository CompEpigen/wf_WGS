#!/usr/bin/env python3


import os
import sys
import numpy as np
import pandas as pd


d_CNA={"sample":[],"chr":[],"start":[],"end":[],"cn":[]}
for x in sorted(os.listdir(sys.argv[1])):
    patient = x.rstrip("_CNVs")
    print(x)
    if os.stat(os.path.join(os.path.join(sys.argv[1],x))).st_size==0:
        continue

    df_CNA = pd.read_csv(os.path.join(os.path.join(sys.argv[1],x)),sep="\t",header=None)
    for i in df_CNA.index:
        chr=df_CNA.iloc[i,0]
        start=df_CNA.iloc[i,1]
        end=df_CNA.iloc[i,2]
        cn=df_CNA.iloc[i,3]
        if end-start>50000:
            d_CNA["sample"].append(patient)
            d_CNA["chr"].append(chr)
            d_CNA["start"].append(start)
            d_CNA["end"].append(end)
            d_CNA["cn"].append(cn)

df_CNA=pd.DataFrame(d_CNA)
df_CNA.to_csv("CNAs.tsv",sep="\t",index=False)
