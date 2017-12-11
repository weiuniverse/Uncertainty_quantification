import math
import tsv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root_path = "../data/poly_mo/"

# preprocess data files

# Quant_bootstraps.tsv :containing the matrix of bootstrap experiments 
# containing the final count for each transcript in each round of bootstrapping 
# with a row be a bootstrap output and columns be list of transcripts. 
quant_bootstraps = tsv.TsvReader(open(root_path+"quant_bootstraps.tsv"))
count = 0
quant_boot = []
for parts in quant_bootstraps:
    quant_boot.append(parts)
df_quant_boot = pd.DataFrame.from_records(quant_boot[1:], columns=quant_boot[0])
df_quant_boot = df_quant_boot.astype('float')
df_quant_boot_mean = df_quant_boot.mean()
df_quant_boot_std = df_quant_boot.std()
id_in_quant_boot = list(df_quant_boot.columns)

# given serial distance
distance

for tid in id_in_quant_boot:
    if tid in distance:
        for i in range(len(df_quant_boot[tid])):
            if -1000 < distance[tid] < 1000:
                df_quant_boot[tid][i] -= distance[tid]