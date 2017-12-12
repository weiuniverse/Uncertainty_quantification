'''
INPUT: read from "quant_bootstraps.tsv", "quant.sf";
OUTPUT: write modified bootstrap matrix to file "features.data"
'''
import tsv
import pandas as pd

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
percent2dot5 = df_qb_sorted.loc[int(sum*0.025)-1]
percent97dot5 = df_qb_sorted.loc[int(sum*0.975)-1]

remove_ids = []
for id in columns:
    if float(percent97dot5[id]) == 0:
        remove_ids.append(id)

new_columns = list(set(columns) - set(remove_ids))
fixed_dfquant_boot = fixed_dfquant_boot[new_columns] 
fixed_dfquant_boot = fixed_dfquant_boot.astype('float')

label = []
for i in range(len(new_columns)):
    label.append(0)
    
nc = [new_columns,label]
nc = list(map(list,zip(*nc)))

new_col = pd.DataFrame.from_records(nc, columns=['Name','useless'])

fixed_dfquant = new_col.merge(df_quant, on='Name')
useless = fixed_dfquant.pop('useless')

df_attributes = fixed_dfquant.copy()
fixed_qb_mean = fixed_dfquant_boot.mean()
fixed_qb_std = fixed_dfquant_boot.std()
fixed_qb_mean = fixed_qb_mean.tolist()
fixed_qb_std = fixed_qb_std.tolist()

df_attributes.insert(5,'quant_boot_mean', fixed_qb_mean)
df_attributes.insert(6,'quant_boot_std', fixed_qb_std)

# write data to file
import pickle
filepath = "features.data"
input_data = df_attributes.as_matrix()
with open(filepath, 'wb+') as output_file:
    pickle.dump(input_data, output_file)

# to read:
# with open(filepath, "rb") as input_file:
#     input_data = pickle.load(input_file)
