#!/usr/bin/env python3


# From a VCF of SVs generated by manta, filter out SVs:
# - whose length is too short
# - which have insufficient Read Pair or Split Reads support
# - (optionally) for which one breakend is in a region present in the Panel of Normal
# - short insertions which likely correspond to retrotransposons (except when they are in the vicinity of a gene). This filtering requires 2 breakpoints.
# (optionally) flags SVs for which one breakend is close to a low mappability region
# (optionally) checks that deletions/duplications are supported by the coverage

import os
import sys
import argparse
import numpy as np
import vcfpy
import pysam
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, help='VCF generated by manta')
parser.add_argument('-o', type = str, help='Output VCF')
parser.add_argument('--pon', type = str, help='PoN containing SVs to filter out.')
parser.add_argument('--minPR', type = int,default=0, help='Minimum number of paired reads supporting the SV')
parser.add_argument('--minSR', type = int,default=0, help='Minimum number of split reads supporting the SV')
parser.add_argument('--minLen', type = int,default=20000, help='Minimum SV length')
parser.add_argument('--bam', type = str, help='BAM file (used to check the breakpoints')
parser.add_argument('--mappability', type = str, help='bed file containing the regions with low mappability')
parser.add_argument('--cnv', type = str, help='cnv file produced by control freec. Less stringent filters are used for SVs close to CNA borders.')
parser.add_argument('--tumorindex', type = int,default=0, help='Index of the tumor sample (0 if only the tumor sample was provided, 1 if a control was used).')
parser.add_argument('--keepCloseSV', type = int,default=0, help='If set to 1, will not filter out two SVs whose breakpoints are close to each other (could be reciprocal translocation)')
parser.add_argument('--filterSmallInsertions', type = int,default=1, help='If set to 0, will not filter out SVs resulting in small insertions.')
parser.add_argument('--puretumor', type = int,default=0, help='If set to 1, will assume that the tumor is pure and will not filter out when no reads support the reference allele..')
parser.add_argument('--log', type = str, help='File where the log will be stored')
args = parser.parse_args()


#sample = "C010-AML-15PB8708"
#args.i = "/home/e840r/Documents/WGS/SVs/manta/filtered/"+sample+".vcf"
#args.o = "/home/e840r/Documents/WGS/SVs/manta/filtered_PoN/"+sample+".vcf"
#args.f= "/home/e840r/Documents/WGS/healthy/SV/PoN_SV.bed"

window_size = 3000 # for filtering out small insertions. 

chromosomes = ["chr"+str(x) for x in range(1,23)] + [str(x) for x in range(1,23)] +["X","Y","chrX","chrY"]
if args.log is not None:
    log_file = open(args.log,"w")
    sys.stdout = log_file


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
positions = set()
reader = vcfpy.Reader.from_path(args.i)
for record in reader:
    chr = str(record.CHROM)
    pos = record.POS
    if not chr in chromosomes: continue
    positions.add((chr,pos))
    if record.INFO["SVTYPE"][:3] in ["DEL","DUP","INS","INV"]:
        chr2 = str(chr)
        pos2 = record.INFO["END"]
    else:
        chr2 = str(record.ALT[0].mate_chrom)
        pos2 = record.ALT[0].mate_pos
    if not chr2 in chromosomes: continue
    positions.add((chr2,pos2))
    breakpoints.add((chr,pos,chr2,pos2))

def sort_chr_pos(t):
    chr = t[0]
    pos = t[1]
    if chr.isdigit(): chr = int(chr)
    else:
        if chr=="X": chr = 23
        else: chr=24
    return chr * 400000000 + pos

positions = sorted(list(positions), key = sort_chr_pos)
breakpoints = sorted(list(breakpoints), key = sort_chr_pos)
def round_position(x):
    return 100*(x//100)

def bps_close(bp1,bp2,margin=200):
    return bp1[0]==bp2[0] and bp1[2]==bp2[2] and abs(bp1[1]-bp2[1])<=margin and abs(bp1[3]-bp2[3]) <=margin


breakpoints_chr_grouped={}
for (chr,pos,chr2,pos2) in breakpoints:
    if (not chr in chromosomes) or (not chr2 in chromosomes): continue
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
        #print(line)
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


def collect_SVs_close(chr,pos,w):
    records= []
    reader_unfiltered = vcfpy.Reader.from_path(args.i)
    for record in reader_unfiltered:
        if chr == record.CHROM and abs(pos-record.POS) < w and pos!=record.POS:
            records.append(record)
    return records


if args.bam is not None:
    samfile = pysam.AlignmentFile(args.bam, "rb")
else:
    print("BAM file was not provided")



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



def low_mapq_in_region(samfile,chr,start,end,mapq_threshold=45):
    count_lowmapq=0
    count_highmapq=0
    for read in samfile.fetch(chr,start,end):
        if read.mapq<mapq_threshold:
            count_lowmapq+=1
        else:
            count_highmapq+=1
    print("low mapq in region " + str(count_lowmapq)+"; high "+str(count_highmapq))
    if pos_close_to_CNA(chr,(start+end)/2):
        return count_lowmapq/(count_lowmapq+count_highmapq) >=0.40
    else:
        return count_lowmapq/(count_lowmapq+count_highmapq) >=0.25

def low_mapq_SV(samfile,chr,pos,chr2,pos2,min_reads,mapq_threshold=45):
    count_supporting_highMAPQ=0
    count_supporting_lowMAPQ=0
    for read in samfile.fetch(chr,pos-200,pos+200):
        if read.next_reference_name==chr2 and abs(read.next_reference_start-pos2)<1000:
            if read.mapq<mapq_threshold:
                count_supporting_lowMAPQ+=1
            else:
                count_supporting_highMAPQ+=1
    print("Read supporting the SV: low mapq " + str(count_supporting_lowMAPQ)+"; high mapq " + str(count_supporting_highMAPQ))
    return count_supporting_highMAPQ<min_reads or count_supporting_highMAPQ*2 < count_supporting_lowMAPQ

def SV_explained_by_alternative_alignments(samfile,chr1,pos1,chr2,pos2):
    # Find reads supporting the alignment
    count_reads_supportingSV=0
    count_reads_supportingSV_withXA=0
    
    for read in samfile.fetch(chr1, pos1-200, pos1+200):
        # Check for split-read (supplementary alignment)
        SA_chr=""
        for t in read.tags:
            if t[0]=="SA":
                for s in t[1][:-1].split(";"):
                    s_split = s.split(",")
                    s_chr = s_split[0]
                    s_start = int(s_split[1])
                    if (s_chr==chr2 and abs(pos2-s_start)<600):
                        SA_chr=s_chr
                        SA_start = s_start
        # Check for read pair supporting the SV
        read_pair_SV = read.next_reference_name==chr2 and abs(read.next_reference_start-pos2)<1000

        if SA_chr!="" or read_pair_SV: # The read supports the SV
            count_reads_supportingSV+=1
            for t in read.tags:
                if t[0]=="XA":
                    count_reads_supportingSV_withXA+=1

    print("Reads supporting SV: " +str(count_reads_supportingSV)+"; with XA: "+str(count_reads_supportingSV_withXA))
    return (count_reads_supportingSV_withXA>0.5*count_reads_supportingSV or (count_reads_supportingSV-count_reads_supportingSV_withXA)<args.minPR)


def filter_min_SR_PR(samfile,chr1,pos1,chr2,pos2):
    # Check that we have enough split reads and paired reads supporting the SV, ignoring those with low mapqual or with alternative alignments
    reads_SR=set()
    reads_PR=set()
    reads_filtered=set()
    mapq1_PR={}
    for read in samfile.fetch(chr1, pos1-200, pos1+200):
        if read.is_supplementary: continue
        has_XA=False
        has_SA=False
        for t in read.tags:
                if t[0]=="XA": has_XA=True
                elif t[0]=="SA": has_SA=True
        if (has_XA and not has_SA) or (read.mapq<=35 and not has_SA) or read.mapq<=10:
            reads_filtered.add(read.query_name)
            continue
        # Check for split-read (supplementary alignment)
        has_correct_SA=False
        for t in read.tags:
            if t[0]=="SA":
                for s in t[1][:-1].split(";"):
                    s_split = s.split(",")
                    s_chr = s_split[0]
                    s_start = int(s_split[1])
                    if (s_chr==chr2 and abs(pos2-s_start)<1000):
                        reads_SR.add(read.query_name)
                        has_correct_SA=True
        # Check for read pair supporting the SV
        if read.next_reference_name==chr2 and abs(read.next_reference_start-pos2)<1000 and ((not has_XA) or has_correct_SA):
            reads_PR.add(read.query_name)
            mapq1_PR[read.query_name]=read.mapq

    mapq2_PR={}
    for read in samfile.fetch(chr2, pos2-200, pos2+200):
        if read.is_supplementary: continue
        has_XA=False
        has_SA=False
        for t in read.tags:
                if t[0]=="XA": has_XA=True
                elif t[0]=="SA": has_SA=True
        if (has_XA and not has_SA) or (read.mapq<=35 and not has_SA) or read.mapq<=10:
            reads_filtered.add(read.query_name)
            continue
        # Check for split-read (supplementary alignment)
        has_correct_SA=False
        for t in read.tags:
            if t[0]=="SA":
                for s in t[1][:-1].split(";"):
                    s_split = s.split(",")
                    s_chr = s_split[0]
                    s_start = int(s_split[1])
                    if (s_chr==chr and abs(pos-s_start)<1000):
                        reads_SR.add(read.query_name)
                        has_correct_SA=True
        # Check for read pair supporting the SV
        if read.next_reference_name==chr and abs(read.next_reference_start-pos)<1000 and ((not has_XA) or has_correct_SA):
            reads_PR.add(read.query_name)
            mapq2_PR[read.query_name]=read.mapq

    for x in reads_filtered:
        if x in reads_PR:
            reads_PR.remove(x)
        if x in reads_SR:
            reads_SR.remove(x)
    
    # Make sure that, on each side, at least one read supporting the SV has a high mapq.
    max_mapq1=0
    max_mapq2=0
    for x in mapq1_PR:
        if not x in reads_filtered: max_mapq1 = max(max_mapq1,mapq1_PR[x])
    for x in mapq2_PR:
        if not x in reads_filtered: max_mapq2 = max(max_mapq2,mapq2_PR[x])
        

    print("PR " +str(len(reads_PR))+"; SR "+str(len(reads_SR)))
    print(reads_PR)
    print(reads_SR)
    print("Max mapq: "+str(max_mapq1)+","+str(max_mapq2))
    
    return len(reads_PR)<args.minPR or len(reads_SR)<args.minSR or max_mapq1<50 or max_mapq2<50

def many_mates_in_region(samfile,chr,pos,stringent):
    # In some regions, there are reads whose mates map to many different regions. These regions are unreliable.
    regions_mates={}
    for read in samfile.fetch(chr, pos-20, pos+20):
        if read.mapq>=45:
            new_region=True
            for (chrom,position) in regions_mates:
                if chrom==read.next_reference_name and abs(position-read.next_reference_start)<50000:
                    new_region=False
                    regions_mates[(chrom,position)]+=1
            if new_region:
                regions_mates[(read.next_reference_name,read.next_reference_start)] =1
    print("Regions of mates: ")
    print(regions_mates)
    count_regions=0
    for x in regions_mates:
        if regions_mates[x]>2: count_regions+=1
    if not stringent: # less stringent filters is a CNA is close to the SV.
        return count_regions>=7 or len(regions_mates)>=20
    else:
        return count_regions>=5 or len(regions_mates)>=12

def reads_go_through_insertion(samfile,chr1,pos1,chr2,pos2):
    """Look for read pairs which go through an insertion"""
    count_insertions=0
    for read in samfile.fetch(chr1,pos1-400,pos1+400):
        if readpair_goes_through_insertion(read,chr1,pos1,chr2,pos2,1000):
            count_insertions+=1
    for read in samfile.fetch(chr2,pos2-400,pos2+400):
        if readpair_goes_through_insertion(read,chr2,pos2,chr1,pos1,1000):
            count_insertions+=1
    if count_insertions>0:
        print("Insertion count: "+ str(count_insertions))
    return count_insertions>0


def readpair_goes_through_insertion(read,chr1,pos1,chr2,pos2,window_size=1000):
    """Look for a read pair where both main alignments are in the same region, but there is a supplementary alignment in between which maps to the other region."""
    has_SA=False
    for t in read.tags:
        if t[0]=="SA":
            for s in t[1][:-1].split(";"):
                s_split = s.split(",")
                s_chr = s_split[0]
                s_start = int(s_split[1])
                if (s_chr==chr1 and abs(pos1-s_start)<window_size) or (s_chr==chr2 and abs(s_start-pos2) < window_size):
                    SA_chr=s_chr
                    SA_start = s_start
                    has_SA=True
    if not has_SA: return False

    # Check that both read pairs map to the same region
    if read.next_reference_name != read.reference_name or abs(read.next_reference_start-read.reference_start)>window_size: return False

    # Check that the supplementary alignment is between the main alignment and the mate
    if read.cigarstring.find("S")>0: charS="S"
    if read.cigarstring.find("H")>0: charS="H"
    if read.is_reverse:
        sup_alignment_middle = read.cigarstring.find(charS) < read.cigarstring.find("M")
    else:
        sup_alignment_middle = read.cigarstring.find(charS) > read.cigarstring.find("M")
    if not sup_alignment_middle: return False

    # Check that the main alignment maps to one of the 2 regions and that the supplementary alignment maps to the other region
    if read.reference_name==chr1 and abs(read.reference_start-pos1)<window_size:
        if SA_chr==chr2 and abs(SA_start-pos2)<window_size:
            return True
    elif read.reference_name==chr2 and abs(read.reference_start-pos2)<window_size:
        if SA_chr==chr1 and abs(SA_start-pos1)<window_size:
            return True

    return False

def reads_bordered_by_two_breakpoints(samfile,chr,pos1,pos2):
    """Look for reads which have soft or hard clipping on both ends, and which start and end at the correct locations."""
    pos1,pos2 = min(pos1,pos2),max(pos1,pos2)
    count_insertions=0
    for read in samfile.fetch(chr,pos1-100,pos2+100):
        if abs(read.reference_start-pos1)<=6 and abs(read.query_alignment_length+read.reference_start-pos2)<=6:
            if (read.cigarstring.endswith("S") or read.cigarstring.endswith("H")):
                if ((read.cigarstring.find("S")>0 and read.cigarstring.find("S")<read.cigarstring.find("M")) or (read.cigarstring.find("H")>0 and read.cigarstring.find("H")<read.cigarstring.find("M")) ):
                    count_insertions+=1
    if count_insertions>0:
        print("Insertion count: "+ str(count_insertions))
    return count_insertions>=4
#########################################################

# Read the CNV file
CNA_breakpoints={}
if args.cnv is not None:
    with open(args.cnv,"r") as infile:
        for line in infile:
            linesplit=line.split("\t")
            chr=linesplit[0]
            start = int(linesplit[1])
            end = int(linesplit[2])
            if end-start>100000:
                if not chr in CNA_breakpoints:
                    CNA_breakpoints[chr]=[]
                CNA_breakpoints[chr].append(start)
                CNA_breakpoints[chr].append(end)

def pos_close_to_CNA(chr,pos):
    if chr in CNA_breakpoints:
        for pos2 in CNA_breakpoints[chr]:
            if abs(pos-pos2)<=30000:
                return True
    return False



    
#########################
## Write the filtered VCF
reader = vcfpy.Reader.from_path(args.i)
header = reader.header
#header.add_filter_line({"ID":"MAPPABILITY","Description":"Regions with low mappability close to the breakpoint."})
writer = vcfpy.Writer.from_path(args.o,reader.header)
for record in reader:
    print("-----------")
    print(record)
    ( (chr,pos,orientation) , (chr2,pos2,orientation2) ) = get_breakpoint_info(record)

   
    # Filter based on minimum number of supporting reads
    if (not "PR" in record.calls[args.tumorindex].data) or (not "SR" in record.calls[args.tumorindex].data):
        filtered_n_reads=True
    else:
        filtered_n_reads = (record.calls[args.tumorindex].data["PR"][1] < args.minPR) or (record.calls[args.tumorindex].data["SR"][1] < args.minSR)
        # Somatic SVs (and heterozygous): expect some reads which do not support the SV.
        print((record.calls[args.tumorindex].data["PR"][1],record.calls[args.tumorindex].data["SR"][1]))
        if args.puretumor==0:
            print("Puretumor")
            filtered_n_reads = filtered_n_reads or (record.calls[args.tumorindex].data["PR"][0] < 2) or (record.calls[args.tumorindex].data["SR"][1] < 2)
    if filtered_n_reads:
        print("SV filtered due to insufficient number of supporting reads")
        continue


    # Filter based on minimum SV length
    filtered_length = (chr==chr2) and ( abs(pos-pos2) < args.minLen)
    if filtered_length:
        print("SV filtered due to short length")
        continue

    # Filter based on Panel of Normal
    if (chr,pos,chr2,pos2) in bp_filtered or (chr2,pos2,chr,pos) in bp_filtered:
            print("SV filtered by Panel of Normal")
            continue

    # Filter based on minimum number of supporting reads
    if filter_min_SR_PR(samfile,chr,pos,chr2,pos2):
        print("Filtered out because there were not enough reads supporting the SV, when removing reads with low mapq.")
        continue

    # Filter out SVs where many of the reads in the region have low MAPQ, since these regions are unreliable.
    
    if args.tumorindex==0: # Do not use this filter if we have a control sample.
        if low_mapq_in_region(samfile,chr,pos-100,pos+100,45) or low_mapq_in_region(samfile,chr2,pos2-100,pos2+100,45):
            print("Filtered out because many of the reads near the breakpoint had a low MAPQ.")
            continue

    # Filter out SVs where the supporting reads have low MAPQ
    #if low_mapq_SV(samfile,chr,pos,chr2,pos2,args.minPR) or low_mapq_SV(samfile,chr2,pos2,chr,pos,args.minPR):
    #    print("Filtered out because many reads supporting the SV had a low MAPQ.")
    #    continue

    # Filter out SVs in regions where there are reads whose mates map to many different regions.
    stringent = not (pos_close_to_CNA(chr,pos) or pos_close_to_CNA(chr2,pos2))
    stringent = stringent and args.tumorindex==0
    if many_mates_in_region(samfile,chr,pos,stringent) or many_mates_in_region(samfile,chr2,pos2,stringent):
        print("Filtered out because reads in the region around the SV had mates mapped to many different regions")
        continue


    # Filter out reads with PolyA insertions, since they probably come from retrotransposons
    if "sequence" in dir(record.ALT[0]) and ("AAAAAAAAA" in record.ALT[0].sequence or "TTTTTTTTT" in record.ALT[0].sequence):
        print("PolyA insertion: probably retrotransposon --> filter out")
        continue 


    # See if a second SV can, combined with the first one, lead to a small insertion, which would then be filtered out.
    filter_out_record = False
    if args.tumorindex==0:
        for r in collect_SVs_close(chr,pos,window_size):
            print(r)
            ( (Bchr,Bpos,Borientation) , (Bchr2,Bpos2,Borientation2) ) = get_breakpoint_info(r)
            if Bchr2 ==chr2 and abs(pos2-Bpos2) < window_size:
                # WARNING: Can have reciprocal translocation with small duplication, which would then look like a small insertion... 
                # In principle, we could use the full read information with the BAM, to check... 
                if pos < Bpos:
                    orientationOutward = (orientation=="+")
                    BorientationOutward = (Borientation=="-")
                else:
                    orientationOutward = (orientation=="-")
                    BorientationOutward = (Borientation=="+")
                if pos2<Bpos2:
                    orientationOutward2 = (orientation2=="+")
                    BorientationOutward2 = (Borientation2=="-")
                else:
                    orientationOutward2 = (orientation2=="-")
                    BorientationOutward2 = (Borientation2=="+")
                
                if abs(pos-Bpos)<120 and abs(pos2-Bpos2)<120: 
                    if reads_go_through_insertion(samfile,chr,pos,chr2,pos2):
                        print("Small insertion, because some read pairs going through the whole insertion were detected")
                        filter_out_record=True
                    elif reads_bordered_by_two_breakpoints(samfile,chr,pos,Bpos) or reads_bordered_by_two_breakpoints(samfile,chr2,pos2,Bpos2):
                        print("Small insertion, because some reads have clipping at both breakends.")
                        filter_out_record=True
                elif reads_bordered_by_two_breakpoints(samfile,chr,pos,Bpos) or reads_bordered_by_two_breakpoints(samfile,chr2,pos2,Bpos2):
                    print("Small insertion, because some reads have clipping at both breakends.")
                    filter_out_record=True
                elif abs(pos-Bpos)<1000 and abs(pos2-Bpos2)<1000 and args.keepCloseSV:
                    print("Close SVs") 
                elif (not orientationOutward) and (not BorientationOutward) and (not orientationOutward2) and (not BorientationOutward2):
                    if chr ==chr2:
                        print("Inversion") # TODO: inverted duplication ??
                    else:
                        print("Reciprocal translocation")
                elif orientationOutward and BorientationOutward and orientationOutward2 and BorientationOutward2 and ((abs(pos-Bpos)<=25 and abs(pos2-Bpos2)>130) or (abs(pos2-Bpos2)<=25 and abs(pos-Bpos)>130) or (not "SR" in r.calls[args.tumorindex].data)):
                    # Insertion with TSD. The duplicated part must be small, but the insertion large (otherwise we would have found a read going through the insertion.)
                    print("Small insertion with TSD ")
                    filter_out_record = True
                    continue
                elif ( orientationOutward and BorientationOutward and (not orientationOutward2) and not (BorientationOutward2) and ((not "SR" in r.calls[args.tumorindex].data) or (abs(pos2-Bpos2)<=25 and abs(pos-Bpos)>130)) ) \
                    or ( (not orientationOutward) and (not BorientationOutward) and orientationOutward2 and BorientationOutward2 and ((not "SR" in r.calls[args.tumorindex].data) or (abs(pos-Bpos)<=25 and abs(pos2-Bpos2)>130))):
                    # Insertion from the first side into the second side, with a deletion at the insertion site. The deletion must be <=25bp, and the insertion larger than 130bp, otherwise we would find reads going through the insertion
                    print("Small insertion, with a small deletion at the inserted site -")
                    filter_out_record = True
                    continue
                elif abs(pos-Bpos)<=5: # The orientation might be wrong...
                    if (not orientationOutward) and (not BorientationOutward) and orientationOutward2 and BorientationOutward2 and (( abs(pos2-Bpos2)>130) or (not "SR" in r.calls[args.tumorindex].data)):
                        # Insertion with TSD. The duplicated part must be small, but the insertion large (otherwise we would have found a read going through the insertion.)
                        print("Small insertion with TSD ")
                        filter_out_record = True
                        continue
                    elif (orientationOutward) and (BorientationOutward) and orientationOutward2 and BorientationOutward2 and ((not "SR" in r.calls[args.tumorindex].data) or  abs(pos2-Bpos2)>130):
                        # Insertion from the first side into the second side, with a deletion at the insertion site. The deletion must be <=25bp, and the insertion larger than 130bp, otherwise we would find reads going through the insertion
                        print("Small insertion, with a small deletion at the inserted site -")
                        filter_out_record = True
                        continue
                elif abs(pos2-Bpos2)<=5: # orientation might be wrong.
                    if orientationOutward and BorientationOutward and (not orientationOutward2) and (not BorientationOutward2) and ( abs(pos-Bpos)>130 or (not "SR" in r.calls[args.tumorindex].data)):
                        # Insertion with TSD. The duplicated part must be small, but the insertion large (otherwise we would have found a read going through the insertion.)
                        print("Small insertion with TSD ")
                        filter_out_record = True
                        continue
                    elif orientationOutward and BorientationOutward and orientationOutward2 and BorientationOutward2 and ((not "SR" in r.calls[args.tumorindex].data)and abs(pos-Bpos)>130):
                        # Insertion from the first side into the second side, with a deletion at the insertion site. The deletion must be <=25bp, and the insertion larger than 130bp, otherwise we would find reads going through the insertion
                        print("Small insertion, with a small deletion at the inserted site -")
                        filter_out_record = True
                        continue
                else:
                    print("Rearrangement is not clear -> keep.")
            print("-")

        if filter_out_record and (args.filterSmallInsertions ==1): continue
    
    if SV_explained_by_alternative_alignments(samfile,chr,pos,chr2,pos2) or SV_explained_by_alternative_alignments(samfile,chr2,pos2,chr,pos):
        print("SV explained by alternative alignment.")
        if stringent: continue
        #continue
    print("Keep SV")
    writer.write_record(record)