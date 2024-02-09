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
from figeno import make_figeno



parser = argparse.ArgumentParser()
parser.add_argument('--freec_cnas', type = str, help='Copy number file')
parser.add_argument('--freec_ratios', type = str, help='Copy number ratios for each bin')
parser.add_argument('--purple_cn', type = str, help='Copy number file from purple')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--format', type = str, default="png", help='Output format: png or svg')
parser.add_argument('--genome', type = str,default="hg19", help='Genome version')
parser.add_argument('--sex', type = str,default="F", help='Sex')
args = parser.parse_args()


chromosomes = [str(x) for x in range(1,23)] + ["X"]
if args.sex=="M": chromosomes.append("Y")
for chr in chromosomes:
    config={"general":{"reference":args.genome,"figure_layout":"horizontal"}}
    config["output"] = {"file":args.o+"chr"+chr+"."+args.format,"dpi":400,"width":180.0}
    config["regions"] = [{"chr":chr}]
    config["tracks"] = [{"type":"SV","height": 10.0,"margin_above": 0.0,"bounding_box": True,"vcf":args.sv},
                            {"type":"copynumber","height": 30.0,"margin_above": 0.0,"bounding_box": True,"freec_ratios":args.freec_ratios,
                             "freec_CNAs":args.freec_cnas,"purple_cn":args.purple_cn,"min_cn":None,"max_cn":None},
                            {"type":"chr axis","height": 10.0,"unit":"Mb","margin_above":0.0}]

    make_figeno(config)




