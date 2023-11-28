#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import vcfpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec




parser = argparse.ArgumentParser()
parser.add_argument('--cna', type = str, help='Copy number file')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('--foldback', type = str, help='tsv file containing foldback inversions (format: chr pos orientation)')
parser.add_argument('--sample', type = str, help='sample name')
parser.add_argument('--sex', type = str,default="F", help='Sex: F or M')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--format', type = str, default="png", help='Output format: png or svg')
parser.add_argument('--chrarms', type = str, help='Length of chromosome arms, in order to plot the centromeres.')
parser.add_argument('--nautosomes', type = int,default=22, help='Number of autosomes (22 for humans, 20 for mice)')
#parser.add_argument('--genome', type = str,default="hg19", help='Genome version (hg19 or hg38)')
parser.add_argument('--roundCN', type = int,default=0, help='If 1, round copy numbers to the nearest integer')
args = parser.parse_args()


sample = args.sample
use_foldbacks = args.foldback is not None

genes_highlighted = ["CDK6","MNX1"]
genes_highlighted = ["EPO","EPOR"]
genes_highlighted=[]

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
    
# Read position of centromeres
if args.chrarms is not None:
    centromeres_pos={}
    with open(args.chrarms,"r") as infile:
        for line in infile:
            linesplit = line.split(" ")
            if linesplit[0][-1]=="p":
                centromeres_pos[linesplit[1]] = int(linesplit[4])

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['font.sans-serif'] = "helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
df_CNA = pd.read_csv(args.cna,sep="\t",dtype={"chromosome":str})
#df_SV = pd.read_csv(args.sv,sep="\t",header=0,dtype={"chr":str})

# Read SV
d_SV={"chr":[],"pos":[],"chr2":[],"pos2":[],"type":[],"color":[]}
reader = vcfpy.Reader.from_path(args.sv)
for record in reader:
    chr = record.CHROM
    pos = record.POS
    chr2 = record.ALT[0].mate_chrom
    pos2 = record.ALT[0].mate_pos
    orientation = record.ALT[0].orientation
    orientation2 = record.ALT[0].mate_orientation
    if chr==chr2 and pos>pos2: continue
    d_SV["chr"].append(chr)
    d_SV["pos"].append(pos)
    d_SV["chr2"].append(chr2)
    d_SV["pos2"].append(pos2)
    if chr!=chr2:
        d_SV["type"].append("TRANS")
        d_SV["color"].append("#27ae60")
    elif orientation=="-" and orientation2=="+":
        d_SV["type"].append("DEL")
        d_SV["color"].append("#4a69bd")
    elif orientation=="+" and orientation2=="-":
        d_SV["type"].append("DUP")
        d_SV["color"].append("#e55039")
    elif orientation=="-" and orientation2=="-":
        d_SV["type"].append("H2H")
        d_SV["color"].append("purple")
    else:
        d_SV["type"].append("T2T")
        d_SV["color"].append("#e58e26")
df_SV = pd.DataFrame(d_SV)

if use_foldbacks:
    df_foldback = pd.read_csv(args.foldback,sep="\t",header=None,names=["chr","pos","orientation"],dtype={"chr":str})


chromosomes = [str(x) for x in range(1,args.nautosomes+1)] + ["X"]
if args.sex=="M":
    chromosomes = chromosomes + ["Y"]
else:
    open(args.o+"_chrY."+args.format,"w").close() # create empty chrY file so that the pipeline knows that all output files are there.

for chr in chromosomes:
    print(chr)
    df_CNA_chr = df_CNA.loc[df_CNA["chromosome"]==chr,:].copy().reset_index()
    df_SV_chr = df_SV.loc[df_SV["chr"]==chr,:]
    df_SV_chr.reset_index(inplace=True)
    if use_foldbacks:
        df_foldback_chr = df_foldback.loc[df_foldback["chr"]==chr,:]
        df_foldback_chr.reset_index(inplace=True,drop="index")

    # Set axes
    fig = plt.figure(figsize=(20,10))
    spec = gridspec.GridSpec(ncols=1, nrows=2,hspace=0.0, height_ratios=[1, 3])
    ax=[]
    ax.append(fig.add_subplot(spec[0]))
    ax.append(fig.add_subplot(spec[1]))

    ymax = np.max([2.6,np.max(df_CNA_chr["copyNumber"])*1.1])
    ymin = np.min([0.9,np.min(df_CNA_chr["copyNumber"])*0.9])
    xmax = np.max(df_CNA_chr["end"])*1.02
    print(xmax)

    if args.chrarms is not None:
        ax[1].vlines(centromeres_pos[chr],ymin,ymax,color="black",linewidth=2)
    

    ax[0].set_xlim(0,xmax)
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



    chr_length = np.max(df_CNA_chr["end"])
    if chr_length>120000000:
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(20000000))
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(2000000))
    else:
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(10000000))
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(1000000))
    if ymax<10:
        ax[1].yaxis.set_major_locator(plt.MultipleLocator(1))
        ax[1].yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    else:
        ax[1].yaxis.set_major_locator(plt.MultipleLocator(5))
        ax[1].yaxis.set_minor_locator(plt.MultipleLocator(1))
    def format_func(value, tick_number):
            return("{:.0f}".format(value/1000000))

    def format_func2(value, tick_number):
        if value%50000000==0:
            return("{:.0f}".format(value/1000000))
        else:
            return ""

    ax[1].ticklabel_format(axis="x", style="sci", scilimits=(6,6))
    ax[1].xaxis.set_major_formatter(plt.FuncFormatter(format_func))

    # CNA plot
    for x in df_CNA_chr.index:
        expected_CN=1 if chr in ["X","Y"] and args.sex=="M" else 2
        if df_CNA_chr.loc[x,"copyNumber"]<0.8*expected_CN: color = "#4a69bd"
        elif df_CNA_chr.loc[x,"copyNumber"]>1.2*expected_CN: color = "#e55039"
        elif df_CNA_chr.loc[x,"bafCount"]>30 and df_CNA_chr.loc[x,"baf"]>0.88: color = "#f6b93b"
        else: color = "black"
        start = df_CNA_chr.loc[x,"start"]
        end = df_CNA_chr.loc[x,"end"]
        if end-start<300000:
            s = 300000-end+start
            start = start-s/2
            end = end+s/2
        linewidth = max(5,40/(ymax-ymin))
        cn = df_CNA_chr.loc[x,"copyNumber"]
        if args.roundCN==1: cn = round(cn)
        ax[1].hlines(cn,start,end,linewidth=linewidth,color=color)
        #ax[1].plot([df_CNA_chr.loc[x,"start"],df_CNA_chr.loc[x,"end"]],[df_CNA_chr.loc[x,"copyNumber"],df_CNA_chr.loc[x,"copyNumber"]],linewidth=5,color=color)

    #ax[1].hlines(1.3,30000000,30001000,color="purple")

    #Highlight some genes
    #for gene_name in genes_highlighted:
    #    gene = data.genes_by_name(gene_name)[0]
    #    if gene.contig!=chr: continue
    #    gene_pos = (gene.start + gene.end) //2
    #    # Find copy number corresponding to the gene (for the y coordinate)
    #    cn_gene=2
    #    for i in range(df_CNA_chr.shape[0]):
    #        if gene_pos >=df_CNA_chr.loc[i,"start"] and gene_pos <=df_CNA_chr.loc[i,"end"]:
    #            cn_gene = df_CNA_chr.loc[i,"copyNumber"] 
    #           if args.roundCN==1: cn_gene = round(cn_gene)
    #    q = ax[1].plot([gene_pos],[cn_gene],marker="x", markersize=20, markeredgecolor="#333333", markerfacecolor="green")
    #    ax[1].text(gene_pos,cn_gene + 0.025*(ymax-ymin)+0.06,gene_name,horizontalalignment="center",style='italic',color="#333333")


    # SV 
    max_SV_size= 10000000 # 1000
    for i in range(df_SV_chr.shape[0]):
        if df_SV_chr.loc[i,"chr"]==df_SV_chr.loc[i,"chr2"] and abs(df_SV_chr.loc[i,"pos"] - df_SV_chr.loc[i,"pos2"]) > max_SV_size:
            max_SV_size = df_SV_chr.loc[i,"pos2"] - df_SV_chr.loc[i,"pos"]
    for i in range(df_SV_chr.shape[0]):
        print(i)
        height = 0.2 + 1.3 * (df_SV_chr.loc[i,"pos2"] - df_SV_chr.loc[i,"pos"]) / max_SV_size
        if df_SV_chr.loc[i,"type"]=="TRANS":
            top = 0.2 + 1.6 * (chr_to_int(df_SV_chr.loc[i,"chr2"]))/27
            ax[0].vlines(x=df_SV_chr.loc[i,"pos"],ymin=0,ymax=top,color="green")
            ax[0].text(x=df_SV_chr.loc[i,"pos"],y=top,s=df_SV_chr.loc[i,"chr2"],horizontalalignment="center",color="black",fontsize=14)
        else:
            add_half_ellipse(ax[0],df_SV_chr.loc[i,"pos"],df_SV_chr.loc[i,"pos2"],height,color=df_SV_chr.loc[i,"color"])
    # Foldback inv
    if use_foldbacks:
        for i in range(df_foldback_chr.shape[0]):
            ax[0].vlines(x=df_foldback_chr.loc[i,"pos"],ymin=0,ymax=1,color="gray")
            if df_foldback_chr.loc[i,"orientation"]=="left":
                ax[0].arrow(x=df_foldback_chr.loc[i,"pos"],y=0.5,dx=-0.008*xmax,dy=0,color="gray",width=0.01,head_width=0.07,head_length=0.005*xmax,length_includes_head=True)
            else:
                ax[0].arrow(x=df_foldback_chr.loc[i,"pos"],y=0.5,dx=+0.008*xmax,dy=0,color="gray",width=0.01,head_width=0.07,head_length=0.005*xmax,length_includes_head=True)

    fig.align_ylabels(ax)
    fig.savefig(args.o+"_chr"+chr+"."+args.format,bbox_inches="tight",pad_inches=0.1,dpi=200)
    plt.cla() 
    plt.clf() 
    plt.close('all')
    #plt.show()






