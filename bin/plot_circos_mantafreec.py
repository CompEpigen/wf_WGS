#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import vcfpy
import pycircos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

parser = argparse.ArgumentParser()
parser.add_argument('--ratios', type = str, help='Copy number file')
parser.add_argument('--cnv', type = str, help='Copy number file')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('--data', type = str, help='Data dir')
parser.add_argument('--sex', type = str,default="XX", help='Sex: XX or XY')
parser.add_argument('-o', type = str, help='Output file')
args = parser.parse_args()

selected_chromosomes = [str(x) for x in range(1,23)] + ["X","Y"]

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



data_dir = args.data

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


color_map={"green":"#27ae60","blue":"#4a69bd","red":"#e55039","black":"black"}
CNV_provided = (args.cnv is not None and os.stat(args.cnv).st_size != 0)
if CNV_provided:
    df_CNV = pd.read_csv(args.cnv,sep="\t",header=None)
    df_CNV.columns=["chr","start","end","cn","type"]
    df_CNV["chr"] = [str(x) for x in df_CNV["chr"]]

ploidy=2
#scatter plot
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open(args.ratios,"r") as f:
    f.readline()
    nskipped=0
    groupsize=0
    last_color="no"
    last_value=-1
    last_name="-1"
    for line in f:
        line  = line.rstrip().split("\t")
        name  = line[0]     
        start = int(line[1])-1
        mid   = start+5000
        #if args.sex =="XY" and name in ["chrX","chrY"]:
        if False:
            value = min(float(line[2]),3.9)
        else:
            value = min(float(line[2])*ploidy,ploidy+1.9)
        if value<0: continue
        values_all.append(value)
        color="black"
        if CNV_provided:
            for i in df_CNV.index:
                if line[0]==df_CNV.loc[i,"chr"] and mid>=df_CNV.loc[i,"start"] and mid<=df_CNV.loc[i,"end"]:
                    if df_CNV.loc[i,"type"]=="gain": color = "#e55039"
                    else: color = "#4a69bd"
        else:
            color="#4a69bd"
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["values"] = []
            arcdata_dict[name]["colors"] = []
        same_as_last = abs(last_value-value) < 0.5 and value > 0.5 and value < 5 and last_name ==last_name and last_color ==color
        if same_as_last and groupsize>3 and nskipped<=20: # skip some points, when a large group is similar...
            nskipped+=1
            groupsize+=1
        else:
            arcdata_dict[name]["positions"].append(mid) 
            arcdata_dict[name]["values"].append(value)
            arcdata_dict[name]["colors"].append(color)
            if same_as_last:
                groupsize+=1
                nskipped=0
            else:
                groupsize=0
            last_name=name
            last_color=color
            last_value = value
    
vmin, vmax = min(values_all), max(values_all) 
vmin2=vmin-0.05*abs(vmin)
vmax2= vmax+0.05*abs(vmax)

for chr in chr_lengths:
    if chr!="chrY":
        for y in [1,2,3]:
            circle.lineplot(chr,data=[y]*20,positions=np.linspace(0,chr_lengths[chr],20),rlim=[vmin2, vmax2], raxis_range=(500,740),linecolor="grey",linewidth=0.5)

for key in arcdata_dict:
    circle.scatterplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                       rlim=[vmin2,vmax2], raxis_range=(500,740), facecolor=arcdata_dict[key]["colors"], spine=True) #scatter plot


#linkplot
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)

d_SV={"chr":[],"pos":[],"chr2":[],"pos2":[],"type":[],"color":[]}
reader = vcfpy.Reader.from_path(args.sv)
for record in reader:
    ( (chr,pos,orientation) , (chr2,pos2,orientation2) ) = get_breakpoint_info(record)
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
    source = (chr1, df_SV.loc[x,"pos"],df_SV.loc[x,"pos"]+1, 500)
    destination =  (chr2, df_SV.loc[x,"pos2"],df_SV.loc[x,"pos2"]+1, 500)
    linewidth=0.8
    if color=="green": linewidth=2.0 # purple
    if chr1 in selected_chromosomes and chr1 in selected_chromosomes:
        circle.chord_plot(source, destination, edgecolor=color,linewidth=linewidth)
circle.figure

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
r1=500+240*(0.87-vmin2)/(vmax2-vmin2)
theta1,r1 = offset(theta1,r1,-0,-25)
theta2= 0
r2=500+240*(1.87-vmin2)/(vmax2-vmin2)
theta2,r2 = offset(theta2,r2,-0,-25)
theta3= 0
r3=500+240*(2.87-vmin2)/(vmax2-vmin2)
theta3,r3 = offset(theta3,r3,0,-25)

theta= 0
r=500+240*0.106
theta,r = offset(theta,r,0,-60)

plt.text(theta1,r1,"1",rotation=0,fontsize=9)
plt.text(theta2,r2,"2",rotation=0,fontsize=9)
plt.text(theta3,r3,"3",rotation=0,fontsize=9)
plt.text(theta,r,"Copy number",rotation=90,fontsize=9)

plt.savefig(args.o,dpi=400,bbox_inches="tight",pad_inches=0)
