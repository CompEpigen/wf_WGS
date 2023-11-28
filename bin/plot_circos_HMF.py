#!/usr/bin/env python3
import argparse
import os
import vcfpy
import pycircos
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['font.sans-serif'] = "helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"

Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(phi, rho)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def offset(theta,r,wx,wy):
    x,y=pol2cart(r,theta)
    return cart2pol(x+wx,y+wy)

parser = argparse.ArgumentParser()
parser.add_argument('--cna', type = str, help='Copy number file')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('--data', type = str, help='Data dir')
parser.add_argument('--sex', type = str,default="XX", help='Sex: XX or XY')
parser.add_argument('-o', type = str, help='Output file')
args = parser.parse_args()
args.data = "/home/e840r/Documents/Scripts/pipelines/WGS_pipeline/data"
data_dir = args.data









selected_chromosomes = [str(x) for x in range(1,23)] + ["X","Y"]

#Set chromosomes
#color_arcs = ["#00FEFD","#FF8A82","#83AEFB","#FCFF44","#4D87FA","#FEE587","#E640F6","#A9FF69","#BBA9A3","#834AF9","#FF5658","#AEBCC3","#83AEFB","#F1FF8A",
#              "#FF80AA","#CE151D","#2F7F1E","#FF9931","#C9FFFE","#CE151D"]
color_arcs=["#98671F","#65661B","#969833","#CE151D","#FF1A25","#FF0BC8","#FFCBCC","#FF9931","#FFCC3A","#FCFF44","#C4FF40","#00FF3B",
            "#2F7F1E","#2800C6","#6A96FA","#98CAFC","#00FEFD","#C9FFFE","#9D00C6","#D232FA","#956DB5","#5D5D5D","#989898","#CBCBCB"]
color_chr={}
for i,chr in enumerate([str(x) for x in range(1,23)]+["X","Y"]):
    color_chr[chr] = color_arcs[i%len(color_arcs)]
circle = Gcircle(figsize=(8,8)) 
chr_lengths={}
df_data_general = pd.read_csv(data_dir+"/circos_chromosome_general.csv",index_col="chr")
df_data_general.index = [x.lstrip("chr") for x in df_data_general.index]

for chr in selected_chromosomes:
    length = df_data_general.loc[chr,"end"]
    chr_lengths[chr] = length
    arc    = Garc(arc_id=chr, size=length, interspace=1, raxis_range=(750,780), labelposition=53, labelsize=18,label_visible=True,facecolor=color_chr[chr],linewidth=0.1) # ,facecolor="white"
    if (args.sex!="XX" or chr!="Y") and chr in selected_chromosomes:
        circle.add_garc(arc) 

circle.set_garcs(00,353) 


#cytoband
import collections
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#D82322"}
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#00000077"}

arcdata_dict = collections.defaultdict(dict)
with open(data_dir+"/circos_chromosome_centromeres.csv") as f: #cytoband
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0].lstrip("chr")
        start = int(line[1])-1 
        width = int(line[2])-(int(line[1])-1) 
        if (args.sex!="XX" or name!="Y") and name in selected_chromosomes:
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["widths"]    = [] 
                arcdata_dict[name]["colors"]    = [] 
            arcdata_dict[name]["positions"].append(start) 
            arcdata_dict[name]["widths"].append(width)
            arcdata_dict[name]["colors"].append(color_dict[line[-1]])

for key in arcdata_dict:
    circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], raxis_range=[750,780], facecolor=arcdata_dict[key]["colors"])    


#blue  #4a69bd

color_map={"green":"#27ae60","blue":"#6AACD2","red":"#e55039","black":"black","purple":"purple"}
color_map_edges={"green":"gpick","blue":"#2E447F","red":"#BD2F19","black":"#292929","purple":"purple"}
#color_map_edges={"green":"gpick","blue":"pink","red":"pink","black":"pink"}

#CN plot
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
df_cna = pd.read_csv(args.cna,sep="\t")
cn_max = 3.8
cn_min = 0.0

for chr in chr_lengths:
    if not chr in selected_chromosomes: continue
    for y in [1,2,3]: #4,5
        circle.lineplot(chr,data=[y]*20,positions=np.linspace(0,chr_lengths[chr],20),rlim=[cn_min, cn_max], raxis_range=(500,740),linecolor="grey",linewidth=0.5)


for x in df_cna.index:
    chr = df_cna.loc[x,"chromosome"]
    if not chr in selected_chromosomes: continue
    cn = df_cna.loc[x,'copyNumber']
    if cn >cn_max: cn = cn_max
    color="black"
    if args.sex=="XX" or (not chr in ['X','Y']):
        if cn<1.5: color = "#4a69bd"
        elif cn>2.4: color = "#e55039"
    else:
        if cn<0.5: color = "#4a69bd"
        elif cn>1.4: color = "#e55039"

    start = max(4000000,df_cna.loc[x,"start"])
    end = min(chr_lengths[chr]-4000000,df_cna.loc[x,'end'])


    circle.lineplot(chr, data=[cn-0.05]*20, positions=np.linspace(start,end,20), 
                    rlim=[cn_min, cn_max], raxis_range=(500,740), linecolor=color, linewidth=3,spine=True)
    circle.lineplot(chr, data=[cn+0.05]*20, positions=np.linspace(start,end,20), 
                    rlim=[cn_min, cn_max], raxis_range=(500,740), linecolor=color, linewidth=3,spine=True)


#linkplot


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

for x in df_SV.index:
    chr1= df_SV.loc[x,"chr"]
    chr2= df_SV.loc[x,"chr2"]
    color = df_SV.loc[x,"color"]
    print(color)
    source = (chr1, df_SV.loc[x,"pos"],df_SV.loc[x,"pos"]+1, 500)
    destination =  (chr2, df_SV.loc[x,"pos2"],df_SV.loc[x,"pos2"]+1, 500)
    print("b")
    linewidth=0.8
    if color=="green": linewidth=2.0 # purple
    if chr1 in selected_chromosomes and chr1 in selected_chromosomes:
        circle.chord_plot(source, destination, edgecolor=color,linewidth=linewidth)
    print("d")

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open(args.sv,"r") as f:
    for line in f:
        line  = line.rstrip().split("\t")
        name1  = line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        color = line[-1].split(",")[0].split("=")[-1]
        source = (name1, start1, end1, 500)
        destination = (name2, start2, end2, 500)
        linewidth=0.8
        if color=="green": linewidth=2.0 # purple
        if name1 in selected_chromosomes and name2 in selected_chromosomes:
            circle.chord_plot(source, destination, edgecolor=color_map[color],linewidth=linewidth)


circle.figure



"""
theta1= 0
r1=500+240*(0.87-vmin2)/(vmax2-vmin2)
theta1,r1 = offset(theta1,r1,-0,-85)
theta2= 0
r2=500+240*(1.87-vmin2)/(vmax2-vmin2)
theta2,r2 = offset(theta2,r2,-0,-85)
theta3= 0
r3=500+240*(2.87-vmin2)/(vmax2-vmin2)
theta3,r3 = offset(theta3,r3,0,-85)

plt.text(theta1,r1,"CN=1",rotation=0,fontsize=8)
plt.text(theta2,r2,"CN=2",rotation=0,fontsize=8)
plt.text(theta3,r3,"CN=3",rotation=0,fontsize=8)

"""

theta1= 0
r1=500+240*(0.87-cn_min)/(cn_max-cn_min)
theta1,r1 = offset(theta1,r1,-0,-25)
theta2= 0
r2=500+240*(1.87-cn_min)/(cn_max-cn_min)
theta2,r2 = offset(theta2,r2,-0,-25)
theta3= 0
r3=500+240*(2.87-cn_min)/(cn_max-cn_min)
theta3,r3 = offset(theta3,r3,0,-25)


theta4= 0
r4=500+240*(3.87-cn_min)/(cn_max-cn_min)
theta4,r4 = offset(theta4,r4,0,-25)
theta5= 0
r5=500+240*(4.87-cn_min)/(cn_max-cn_min)
theta5,r5 = offset(theta5,r5,0,-25)

theta= 0
r=500+240*0.106
theta,r = offset(theta,r,0,-60)

plt.text(theta1,r1,"1",rotation=0,fontsize=9)
plt.text(theta2,r2,"2",rotation=0,fontsize=9)
plt.text(theta3,r3,"3",rotation=0,fontsize=9)
#plt.text(theta4,r4,"4",rotation=0,fontsize=9)
#plt.text(theta5,r5,"5",rotation=0,fontsize=9)
plt.text(theta,r,"Copy number",rotation=90,fontsize=9)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
color_map={"green":"#27ae60","blue":"#6AACD2","red":"#e55039","black":"black"}
p= [Line2D([0], [0], color=color_map["green"],label="Translocation",lw=2),
    Line2D([0], [0], color=color_map["blue"],label="Deletion",lw=2),
    Line2D([0], [0], color=color_map["red"],label="Duplication",lw=2),
    Line2D([0], [0], color=color_map["black"],label="Inversion",lw=2)
]
#plt.legend(handles=p,fontsize=12)
#plt.title(sample)


plt.savefig(args.o,dpi=600,bbox_inches="tight",pad_inches=0,transparent=True) #transparent=True
