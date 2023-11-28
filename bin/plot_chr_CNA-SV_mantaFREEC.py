#!/usr/bin/env python3


import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
import vcfpy



parser = argparse.ArgumentParser()
parser.add_argument('--cnv', type = str, help='Copy number file')
parser.add_argument('--ratios', type = str, help='Copy number ratios for each bin')
parser.add_argument('--sv', type = str, help='SV vcf')
parser.add_argument('--foldback', type = str, help='tsv file containing foldback inversions (format: chr pos orientation)')
parser.add_argument('--sample', type = str, help='sample name')
parser.add_argument('--sex', type = str,default="F", help='Sex: F or M')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--format', type = str, help='Output format: png or svg')
parser.add_argument('--nautosomes', type = int,default=22, help='Number of autosomes (22 for humans, 20 for mice)')
args = parser.parse_args()


sample = args.sample
use_foldbacks = args.foldback is not None

def add_half_ellipse(ax,x1,x2,ymax,color):
    e = patches.Ellipse(xy=((x1+x2)/2,0),width=abs(x2-x1),height=ymax*2,fill=False,edgecolor=color)
    ax.add_patch(e)

def chr_to_int(chr):
    chr=str(chr)
    if chr.isdigit():
        return int(chr)
    elif chr=="X":
        return 23
    else:
        return 24

matplotlib.rcParams.update({'font.size': 20})
#matplotlib.rcParams['font.sans-serif'] = "helvetica"
#matplotlib.rcParams['font.family'] = "sans-serif"


CNVs=[]
with open(args.cnv,"r") as infile:
    for line in infile:
        linesplit = line.rstrip("\n").split("\t")
        CNVs.append((linesplit[0],int(linesplit[1])/ 1000000,int(linesplit[2])/ 1000000,linesplit[4]))

def bin2CNV(CNVs,chr,start):
    result="normal"
    for CNV in CNVs:
        if str(chr)==str(CNV[0]) and start>=CNV[1] and start<=CNV[2]:
            result = CNV[3]
    return result

df_ratios = pd.read_csv(args.ratios,sep="\t",dtype={"Chromosome":str,"Start":float})
df_ratios = df_ratios.loc[df_ratios["Ratio"]>=0,:]
#binsize = df_ratios.loc[1,"Start"] - df_ratios.loc[0,"Start"]

if use_foldbacks:
    df_foldback = pd.read_csv(args.foldback,sep="\t",header=None,names=["chr","pos","orientation"],dtype={"chr":str})


chromosomes = [str(x) for x in range(1,args.nautosomes+1)] + ["X"]
if "Y" in args.sex:
    chromosomes = chromosomes + ["Y"]
else:
    open(args.o+"_chrY."+args.format,"w").close() # create empty chrY file so that the pipeline knows that all output files are there.

ploidy=2
for chr in chromosomes:
    print("chr"+str(chr))
    df_ratios_chr = df_ratios.loc[df_ratios["Chromosome"]==chr,:].copy(deep=True)
    df_ratios_chr["Start"] = df_ratios_chr["Start"] / 1000000
    if use_foldbacks:
        df_foldback_chr = df_foldback.loc[df_foldback["chr"]==chr,:]
        df_foldback_chr.reset_index(inplace=True,drop="index")

    # Set axes
    fig = plt.figure(figsize=(20,10))
    spec = gridspec.GridSpec(ncols=1, nrows=2,hspace=0.0, height_ratios=[1, 3])
    ax=[]
    ax.append(fig.add_subplot(spec[0]))
    ax.append(fig.add_subplot(spec[1]))

    ymax = min(np.max([ploidy+1.2,np.max(df_ratios_chr["Ratio"])*ploidy*1.1]),50)
    ymin = max(0,np.min([1.2,np.min(df_ratios_chr["Ratio"])*ploidy*0.9]))
    xmax = np.max(df_ratios_chr["Start"])*1.02
    if xmax!=xmax: xmax=10
    if ymax!=ymax: ymax=3.5
    if ymin!=ymin: ymin=0
    print((xmax,ymin,ymax))
    

    ax[0].set_xlim(0,max(1,xmax))
    ax[0].set_ylim(0,2)
    ax[0].tick_params(axis='y', which='both', left=False,right=False,labelleft=False)
    ax[0].tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)
    ax[0].set_ylabel("Breakpoints")
    ax[0].set_title(sample + " (chr"+str(chr)+")")


    ax[1].set_ylim(ymin,ymax)
    ax[1].set_xlim(0,xmax)
    ax[1].grid(visible=True,which="major",axis="x",color="darkgray")
    ax[1].grid(visible=True,which="minor",axis="x",color="lightgray",linestyle="dashed")
    ax[1].grid(visible=True,which="major",axis="y",color="lightgray",linestyle="dashed")
    ax[1].set_axisbelow(True)
    ax[1].set_xlabel("Position on chr"+str(chr) + " (Mb)")
    ax[1].set_ylabel("Copy number")



    chr_length = np.max(df_ratios_chr["Start"])
    if chr_length>120:
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(20))
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(2))
    else:
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(10))
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(1))
    if ymax<10:
        ax[1].yaxis.set_major_locator(plt.MultipleLocator(1))
        ax[1].yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    else:
        ax[1].yaxis.set_major_locator(plt.MultipleLocator(5))
        ax[1].yaxis.set_minor_locator(plt.MultipleLocator(1))
    def format_func(value, tick_number):
            return("{:.0f}".format(value))

    def format_func2(value, tick_number):
        if value%50==0:
            return("{:.0f}".format(value))
        else:
            return ""

    ax[1].ticklabel_format(axis="x", style="sci", scilimits=(6,6))
    ax[1].xaxis.set_major_formatter(plt.FuncFormatter(format_func))

    # CNA plot
    colors=[]
    for i in df_ratios_chr.index:
        cn_state = bin2CNV(CNVs,df_ratios_chr.loc[i,"Chromosome"],df_ratios_chr.loc[i,"Start"])
        if cn_state=="gain":
            colors.append("#e55039")
        elif cn_state=="loss":
            colors.append("#4a69bd")
        else:
            colors.append("black")
    ax[1].scatter(df_ratios_chr["Start"],ploidy*df_ratios_chr["Ratio"],c=colors,s=2.0,marker="o")


    # SV 
    max_SV_size=1000/1000000

    SVs=[]
    reader = vcfpy.Reader.from_path(args.sv)
    for record in reader:
        chr1 = record.CHROM
        if chr1!=chr: continue
        pos1 = record.POS /1000000
        if record.INFO["SVTYPE"] in ["BND","TRA"]:
            chr2 = record.ALT[0].mate_chrom
            pos2 = record.ALT[0].mate_pos /1000000
        else:
            chr2 = chr1
            pos2 = record.INFO["END"]  /1000000
        
        if record.INFO["SVTYPE"] == "DEL":
            color="#4a69bd"
        elif record.INFO["SVTYPE"] == "DUP":
            color = "#e55039"
        elif record.INFO["SVTYPE"] =="INV":
            color = "purple"
        else:  #BND
            if record.CHROM != record.ALT[0].mate_chrom:
                color = "green"
            else:
                if record.POS > record.ALT[0].mate_pos:
                    continue
                if record.ALT[0].orientation =="-" and record.ALT[0].mate_orientation == "+":
                    color = "#4a69bd"
                elif record.ALT[0].orientation =="+" and record.ALT[0].mate_orientation == "-":
                    color = "#e55039"
                else:
                    color = "purple"
        if chr1==chr2 and abs(pos2-pos1)>max_SV_size: max_SV_size =  abs(pos2-pos1)
        SVs.append((chr1,pos1,chr2,pos2,color))

    for (chr1,pos1,chr2,pos2,color) in SVs:
        height = 0.2 + 1.3 * abs(pos2-pos1) / max_SV_size
        if color=="green":
            top = 0.2 + 1.6 * (chr_to_int(chr2))/27
            ax[0].vlines(x=pos1,ymin=0,ymax=top,color="green")
            ax[0].text(x=pos1,y=top,s=chr2,horizontalalignment="center",color="black",fontsize=14)
        else:
            add_half_ellipse(ax[0],pos1,pos2,height,color="black")

    # Foldback inv
    if use_foldbacks:
        for i in range(df_foldback_chr.shape[0]):
            ax[0].vlines(x=df_foldback_chr.loc[i,"pos"],ymin=0,ymax=1,color="purple")
            if df_foldback_chr.loc[i,"orientation"]=="left":
                ax[0].arrow(x=df_foldback_chr.loc[i,"pos"],y=0.5,dx=-0.008*xmax,dy=0,color="purple",width=0.01,head_width=0.07,head_length=0.005*xmax,length_includes_head=True)
            else:
                ax[0].arrow(x=df_foldback_chr.loc[i,"pos"],y=0.5,dx=+0.008*xmax,dy=0,color="purple",width=0.01,head_width=0.07,head_length=0.005*xmax,length_includes_head=True)

    fig.align_ylabels(ax)
    fig.savefig(args.o+"_chr"+chr+"."+args.format,bbox_inches="tight",pad_inches=0.1,dpi=200)
    plt.cla() 
    plt.clf() 
    plt.close('all')
    #plt.show()






