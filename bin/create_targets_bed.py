#!/usr/bin/env python3

# Create a bed file covering genes, for a specific chromosome. This is used to only look for SNVs in these regions.
# Arguments: gtf, chr, output file

import sys
import gzip

gtf_path = sys.argv[1]
chr = sys.argv[2]
output = sys.argv[3]

# Find genes coordinates in the gtf file
l=[]
if gtf_path.endswith(".gz") : infile = gzip.open(gtf_path,"rt")
else: infile =open(gtf_path,"r")
for line in infile:
    if line.startswith("#"): continue
    linesplit = line.split("\t")
    if linesplit[0].lstrip("chr") != chr: continue
    if linesplit[2]!="gene": continue
    l.append([int(linesplit[3]),int(linesplit[4])])
infile.close()

# Sort coordinates
l = sorted(l)

# Merge overlapping intervals
i=1
while i < len(l):
    prev_start,prev_end = l[i-1]
    if l[i][0] < l[i-1][1]:
        l[i-1] = [l[i-1][0],l[i][1]]
        l.remove(l[i])
    else:
        i+=1

# Write targets to the output file
with open(output,"w") as outfile:
    for x in l:
        tmp = outfile.write("\t".join([chr,str(x[0]),str(x[1])])+"\n")