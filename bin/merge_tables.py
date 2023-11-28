#!/usr/bin/env python3

import pandas as pd
import sys

df_merged=pd.DataFrame()
for f in sys.argv[2:]:
    df = pd.read_csv(f,sep="\t")
    df_merged = pd.concat([df_merged,df])

df_merged.to_csv(sys.argv[1],sep="\t",index=False)
