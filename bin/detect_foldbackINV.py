#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import gzip
import vcfpy

parser = argparse.ArgumentParser()
parser.add_argument('--sv', type = str, help='unfiltered VCF of structural variants')
parser.add_argument('--cn', type = str, help='CN segments')
parser.add_argument('--pon', type = str, help='Panel of normal')
parser.add_argument('-o', type = str, help='Output tsv file containing the foldback inversions')
args = parser.parse_args()

def get_breakpoint_info(record):
    chr = record.CHROM
    pos = record.POS
    chr2 = record.ALT[0].mate_chrom
    pos2 = record.ALT[0].mate_pos
    orientation = record.ALT[0].orientation
    orientation2 = record.ALT[0].mate_orientation
    return ( (chr,pos,orientation) , (chr2,pos2,orientation2) )


def locus_in_bedRegion(file,chr,pos):
    """
    Parse the bed files of filtered regions to see if one position is present in one of these regions.
    Start from the last position, to avoid iterating several times through the file (the files are assumed to be sorted)
    """
    current_chr = "0"
    start = 0
    end=0
    line = True
    #print("------ Looking for "+chr+"_"+str(pos))
  
    last_pos_file = file.tell()
    pos_file = file.tell()
    while line and current_chr!=chr: # get to the right chromosome
        #print(line)
        last_pos_file = file.tell()
        pos_file = file.tell()
        line = file.readline()
        if line:
            linesplit = line.split("\t")
            current_chr = linesplit[0]
            start = int(linesplit[1])
            end = int(linesplit[2])
  
    while line and current_chr==chr and end<pos:
        last_pos_file = pos_file
        pos_file = file.tell()
        line = file.readline()
        if line:
            linesplit = line.split("\t")
            current_chr = linesplit[0]
            start = int(linesplit[1])
            end = int(linesplit[2])
    file.seek(last_pos_file) # go back one line 

    if pos <= end and pos >= start:
        return True
    else:
        return False


def get_CNV_breakpoints(filename):
    #only include breakpoints of gains ? 
    df_CNA = pd.read_csv(filename,sep="\t")
    breakpoints = []
    for x in df_CNA.index:
        if df_CNA.loc[x,"copyNumber"]>2.7:
            if x>0 and abs(df_CNA.loc[x,"copyNumber"]-df_CNA.loc[x-1,"copyNumber"])>0.9:
                breakpoints.append((df_CNA.loc[x,"chromosome"],df_CNA.loc[x,"start"]))
            if x+1<df_CNA.shape[0] and abs(df_CNA.loc[x,"copyNumber"]-df_CNA.loc[x+1,"copyNumber"])>0.9:
                breakpoints.append((df_CNA.loc[x,"chromosome"],df_CNA.loc[x,"end"]))
    return breakpoints

def locus_close_to_CNVbreakpoint(chr,pos,CNVbreakpoints,w=10000):  #20000
    for (chr2,pos2) in CNVbreakpoints:
        if chr==chr2 and abs(pos-pos2)<=w:
            return True
    return False


def detect_foldbackINV(filename,CNVfile,bp_filtered):
    # Based on criteria from https://www.nature.com/articles/s41467-019-13824-9 and https://www.sciencedirect.com/science/article/pii/S000292971500508X?via%3Dihub
    # Inversion where the 2 breakpoints are within 20kb, and no reciprocal partner
    # Also copy number change near the breakpoint
    breakpoints = get_CNV_breakpoints(CNVfile)

    reader = vcfpy.Reader.from_path(os.path.join(filename))
    inversions=[]
    for record in reader:
        ( (chr,pos,orientation) , (chr2,pos2,orientation2) ) = get_breakpoint_info(record)
        if chr==chr2 and orientation == orientation2 and pos<pos2:
            if "RP" in record.calls[0].data: 
                RP = record.calls[0].data["RP"]
            else: 
                RP=0

            if "SR" in record.calls[0].data: 
                SR = record.calls[0].data["SR"]
            else: 
                SR=0
            inversions.append((chr,pos,pos2,orientation,RP,SR))

    foldback_inversions=[]
    for (chr,pos,pos2,orientation,RP,SR) in inversions:
        if RP >=4 and SR >=4 and (not (chr,pos,chr2,pos2) in bp_filtered) and (not (chr2,pos2,chr,pos) in bp_filtered):
            if abs(pos-pos2)<20000: #TODO: papers use 20000 but this seems very large....
                # Check if there is a reciprocal partner
                reciprocal_partner=False
                for (Bchr,Bpos,Bpos2,Borientation,BRP,BSR) in inversions:
                    if chr==Bchr and pos==Bpos and pos2==Bpos2 and orientation==Borientation: continue
                    if chr==Bchr and abs(pos-Bpos)<3000 and abs(pos2-Bpos2)<3000 and orientation!=Borientation:
                        reciprocal_partner=True

                if not reciprocal_partner:
                    # Only consider foldback inversions which are close to a copy number change
                    CNV_close = locus_close_to_CNVbreakpoint(chr,pos,breakpoints) or locus_close_to_CNVbreakpoint(chr,pos2,breakpoints)
                    if CNV_close:
                        if orientation=="-":
                            orientation_renamed="left"
                        else:
                            orientation_renamed="right"
                        foldback_inversions.append((chr,pos,pos2,orientation_renamed))

    return(foldback_inversions)


####################################
# Panel of normal
chromosomes = [str(x) for x in range(1,23)] + ["X"]

def bp_in_pon(file,chr1,pos1,chr2,pos2,margin=5):
    """
    Parse the bed files of filtered regions to see if one position is present in one of these regions.
    Start from the last position, to avoid iterating several times through the file (the files are assumed to be sorted)
    """
    current_chr = "0"
    start = 0
    end=0
    line = True
    #print("------ Looking for "+chr+"_"+str(pos))
  
    last_pos_file = file.tell()
    pos_file = file.tell()
    while line and current_chr!=chr1: # get to the right chromosome
        #print(line)
        last_pos_file = file.tell()
        pos_file = file.tell()
        line = file.readline()
        if line:
            line = line.decode("ascii")
            linesplit = line.split("\t")
            current_chr = linesplit[0]
            start = int(linesplit[1])
            end = int(linesplit[2])
  
    while line and current_chr==chr1 and end+margin<pos1:
        last_pos_file = pos_file
        pos_file = file.tell()
        line = file.readline()
        if line:
            line = line.decode("ascii")
            linesplit = line.split("\t")
            current_chr = linesplit[0]
            start = int(linesplit[1])
            end = int(linesplit[2])
            current_chr2 =  linesplit[3]
            start2 =  int(linesplit[4])
            end2 =  int(linesplit[5])
    file.seek(last_pos_file) # go back one line 

    if pos1 <= end+margin and pos1+margin >= start and current_chr2 ==chr2 and pos2 <= end2+margin and pos2+margin >= start2:
        return True
    else:
        return False



# Get all of the positions of the breakends
breakpoints=set()
reader = vcfpy.Reader.from_path(args.sv)
for record in reader:
    chr = str(record.CHROM)
    pos = record.POS
    if not chr in chromosomes: continue
    chr2 = str(record.ALT[0].mate_chrom)
    pos2 = record.ALT[0].mate_pos
    if not chr2 in chromosomes: continue
    breakpoints.add((chr,pos,chr2,pos2))

def sort_chr_pos(t):
    chr = t[0]
    pos = t[1]
    if chr.isdigit(): chr = int(chr)
    else:
        if chr=="X": chr = 23
        else: chr=24
    return chr * 400000000 + pos

breakpoints = sorted(list(breakpoints), key = sort_chr_pos)
def round_position(x):
    return 100*(x//100)
def bps_close(bp1,bp2,margin=200):
    return bp1[0]==bp2[0] and bp1[2]==bp2[2] and abs(bp1[1]-bp2[1])<=margin and abs(bp1[3]-bp2[3]) <=margin

breakpoints_chr_grouped={}
for (chr,pos,chr2,pos2) in breakpoints:
    pos_rounded = round_position(pos)
    if (not chr in breakpoints_chr_grouped):
        breakpoints_chr_grouped[chr] = {}
    if (not pos_rounded in breakpoints_chr_grouped[chr]):
        breakpoints_chr_grouped[chr][pos_rounded] = []
    breakpoints_chr_grouped[chr][pos_rounded].append((chr,pos,chr2,pos2))



# Go through the PoN file to see which of the breakpoints are filtered out
bp_filtered = set()
if args.pon is not None:
    print("PoN was provided")
    if args.pon.endswith("gz"):
        file_filter = gzip.open(args.pon,"r")
    else:
        file_filter = open(args.pon,"r")
    line = True
    while line: 
        line = file_filter.readline()
        if line:
            line = line.decode("ascii")
            linesplit = line.split("\t")
            chr = linesplit[0]
            start = int(linesplit[1])
            end = int(linesplit[2])
            chr2 =  linesplit[3]
            start2 =  int(linesplit[4])
            end2 =  int(linesplit[5])
            if (not chr in chromosomes) or (not chr2 in chromosomes): continue 
            filtered=False
            positions_to_test=[start]
            if round_position(start-30)!=round_position(start): positions_to_test.append(start-30)
            if round_position(start+30)!=round_position(start): positions_to_test.append(start+30)
            for p in positions_to_test:
                pos_rounded = round_position(p)
                if chr in breakpoints_chr_grouped and pos_rounded in breakpoints_chr_grouped[chr]:
                    for bp in breakpoints_chr_grouped[chr][pos_rounded]:
                        if bps_close((chr,start,chr2,start2),bp):
                            bp_filtered.add(bp)
            
            positions_to_test=[start2]
            if round_position(start2-30)!=round_position(start2): positions_to_test.append(start2-30)
            if round_position(start2+30)!=round_position(start2): positions_to_test.append(start2+30)
            for p in positions_to_test:
                pos_rounded = round_position(p)
                if chr2 in breakpoints_chr_grouped and pos_rounded in breakpoints_chr_grouped[chr2]:
                    for bp in breakpoints_chr_grouped[chr2][pos_rounded]:
                        if bps_close((chr2,start2,chr,start),bp):
                            bp_filtered.add(bp)





foldback_inversions = detect_foldbackINV(args.sv,args.cn,bp_filtered)

with open(args.o,"w") as outfile:
    for inv in foldback_inversions:
        tmp = outfile.write(inv[0]+"\t"+str(inv[1])+"\t"+inv[3]+"\n")


