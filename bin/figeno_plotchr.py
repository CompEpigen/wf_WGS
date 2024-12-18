#!/usr/bin/env python3

import argparse
from figeno import figeno_make



parser = argparse.ArgumentParser()
parser.add_argument('--freec_cnas', type = str, help='Copy number file')
parser.add_argument('--freec_ratios', type = str, help='Copy number ratios for each bin')
parser.add_argument('--purple_cn', type = str, help='Copy number file from purple')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--format', type = str, default="png", help='Output format: png or svg')
parser.add_argument('--genome', type = str,default="hg19", help='Genome version')
parser.add_argument('--sex', type = str,default="F", help='Sex')
parser.add_argument('--ploidy', type = float,default=2.0, help='Ploidy')
args = parser.parse_args()


chromosomes = [str(x) for x in range(1,23)] + ["X"]
if args.sex=="M": chromosomes.append("Y")
for chr in chromosomes:
    config={"general":{"reference":args.genome,"layout":"horizontal"}}
    config["output"] = {"file":args.o+"_chr"+chr+"."+args.format,"dpi":400,"width":180.0}
    config["regions"] = [{"chr":chr}]
    config["tracks"] = [{"type":"sv","height": 10.0,"margin_above": 0.0,"bounding_box": True,"file":args.sv},
                            {"type":"copynumber","height": 30.0,"margin_above": 0.0,"bounding_box": True,"freec_ratios":args.freec_ratios,
                             "freec_CNAs":args.freec_cnas,"purple_cn":args.purple_cn,"ploidy":args.ploidy,"min_cn":None,"max_cn":None},
                            {"type":"chr_axis","height": 10.0,"unit":"Mb","margin_above":0.0}]

    figeno_make(config)




