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

# Quant.sf :estimated quantifications for each transcript
quant_file = open(root_path+"quant.sf")
lines = quant_file.readlines()
quant_file.close()
count = 0
quant = []
for line in lines:
    line = line[:-1]
    l = line.split('\t')
    quant.append(l)

df_quant = pd.DataFrame.from_records(quant[1:], columns=quant[0])
df_quant.Name = df_quant.Name.astype(str)
df_quant.Length = df_quant.Length.astype(int)
df_quant.EffectiveLength = df_quant.EffectiveLength.astype(float)
df_quant.TPM = df_quant.TPM.astype(float)
df_quant.NumReads = df_quant.NumReads.astype(float)

# eq_classes

# we now have:
# df_quant_boot_mean
# df_quant_boot_std
# df_quant.Name
# df_quant.Length
# df_quant.EffectiveLength
# df_quant.TPM
# df_quant.NumReads

fixed_dfquant_boot = df_quant_boot.copy()
columns = fixed_dfquant_boot.columns

sort_qb = []
for id in columns:
    try:
        listed = list(df_quant_boot[id])        
    except KeyError:
        pass
    else:
        listed.sort()
        sort_qb.append(listed)
sort_qb = list(map(list,zip(*sort_qb)))

df_qb_sorted = pd.DataFrame.from_records(sort_qb, columns=columns)

sum = len(sort_qb)
print(int(sum*0.025))
print(int(sum*0.975))
percent2dot5 = df_qb_sorted.loc[int(sum*0.025)-1]
percent97dot5 = df_qb_sorted.loc[int(sum*0.975)-1]

remove_ids = []
for id in columns:
    if float(percent97dot5[id]) == 0:
        remove_ids.append(id)

new_columns = set(columns) - set(remove_ids)
new_columns = list(new_columns)
fixed_dfquant_boot = fixed_dfquant_boot[new_columns] 