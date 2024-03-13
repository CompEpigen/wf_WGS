# Script for creating a gene x sample TPM.tsv file, from a wf_WGS samplesheet with a bam_RNA entry.
# This is intended to be used on the DKFZ cluster.

import pandas as pd
import os
import sys

samplesheet=sys.argv[1]
df_samples = pd.read_csv(samplesheet,sep=",")
firstSample=True
dic={}
for x in df_samples.index:
    sample = df_samples.loc[x,"sample"]
    bam_rna = df_samples.loc[x,"bam_RNA"]
    d = os.path.join(os.path.dirname(bam_rna),"featureCounts")
    for f in os.listdir(d):
        if f.endswith("fpkm_tpm.featureCounts.tsv"):
            df = pd.read_csv(os.path.join(d,f),sep="\t")
            if firstSample:
                gene_ids = [x.split(".")[0] for x in df["gene_id"]]
                gene_names = list(df["name"])
                firstSample=False
            dic[sample] = list(df["TPM"])


df = pd.DataFrame(dic)
cols = df.columns
df["gene_id"] = gene_ids
df["gene_name"] = gene_names
cols = ["gene_id","gene_name"] + list(cols)
df = df[cols]
df.to_csv("TPM.tsv",sep="\t",index=False)