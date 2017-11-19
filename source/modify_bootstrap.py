
# coding: utf-8

# # Uncertainty Quantification

# ## Overview
#     we will analyze why we cannot get the right count for some transcripts using the output of salmon.

# ## Analyze tools
#     we will mainly use dataframe of pandas to analyze the data.


import tsv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
# %matplotlib inline


# # data Preprocess

# ## Poly Truth
#     Read the file poly_truth.tsv
#     Poly_truth.tsv: true counts for each transcript

# Poly_truth.tsv: true counts for each transcript
def read_poly_truth(root_path):
    # read the file
    poly_truth = open(root_path+"poly_truth.tsv")
    lines = poly_truth.readlines()
    poly_truth.close()
    # split every line
    count = 0
    poly_truth = []
    for line in lines:
        line = line[:-1]
        l = line.split('\t')
        poly_truth.append(l)

    df_poly_truth = pd.DataFrame.from_records(poly_truth[1:], columns=poly_truth[0])
    df_poly_truth['transcript_id']=df_poly_truth['transcript_id'].astype(str)
    df_poly_truth['count']=df_poly_truth['count'].astype(int)
    truth_id = df_poly_truth.transcript_id

    return df_poly_truth,truth_id

# ## Quant_bootstraps
#     Read the file quant_bootstraps.tsv
#     Quant_bootstraps.tsv :containing the matrix of bootstrap experiments containing the final count for each transcript in each round of bootstrapping with a row be a bootstrap output and columns be list of transcripts.

# Quant_bootstraps.tsv :containing the matrix of bootstrap experiments
# containing the final count for each transcript in each round of bootstrapping
# with a row be a bootstrap output and columns be list of transcripts.
def read_bootstraps(root_path):
    quant_bootstraps = tsv.TsvReader(open(root_path+"quant_bootstraps.tsv"))
    count = 0
    quant_boot = []
    for parts in quant_bootstraps:
        quant_boot.append(parts)

    df_quant_boot = pd.DataFrame.from_records(quant_boot[1:], columns=quant_boot[0])
    id_qb = list(df_quant_boot.columns)
    return df_quant_boot,id_qb

def sort_bootstrap(df_quant_boot,truth_id):
    sort_qb = []
    use_id = []
    for id in truth_id:
        try:
            listed = list(df_quant_boot[id])
        except KeyError:
            pass
        else:
            use_id.append(id)
            listed.sort()
            sort_qb.append(listed)
    # ### reverse sort_qb
    sort_qb = list(map(list,zip(*sort_qb)))
    # ### transfer to dataframe
    df_qb_sorted = pd.DataFrame.from_records(sort_qb, columns=use_id)

    return df_qb_sorted,sort_qb,use_id





# ## Eq_classes.txt
#     Eq_classes.txt: list of equivalence classes and their information
def read_eq_classes(root_path):
    # read the file
    file = open(root_path+"eq_classes.txt")
    lines = file.readlines()
    file.close()

    # clean the useless data
    counts_transcript = lines[0]
    counts_transcript = int(counts_transcript[:-1])
    counts_class = lines[1]
    counts_class = int(counts_class[:-1])
    lines = lines[2:]

    # ### get the sequence and order of id
    seq_id=[]
    for id in lines[0:counts_transcript]:
        seq_id.append(id[:-1])

    lines = lines[counts_transcript:]

    list_class = []
    for line in lines:
        li = line[:-1].split('\t')
        list_class.append(li)

    list_class = list(list(int(a) for a in b) for b in list_class)
    return list_class,seq_id




# ## Filter the False

# ### find the value of 2.5% & 97.5% and the false transcript is out of this range
def fileter_the_false(df_poly_truth,df_qb_sorted,sort_qb,id_qb,use_id,truth_id,alpha=0.475):
    df_poly_truth = df_poly_truth.set_index(['transcript_id'])
    sum = len(sort_qb)
    low_percent = 0.5-alpha
    high_percent = 0.5+alpha
    low = df_qb_sorted.loc[int(sum*low_percent)-1]
    high = df_qb_sorted.loc[int(sum*high_percent)-1]
    # ## divide the transcript_id into two group
    true_id = []
    false_id = []
    for id in use_id:
        # print(low[id])
        down = float(low[id])
        up = float(high[id])
        true_count = df_poly_truth.loc[id]
        true_count = float(true_count)
        if true_count>down and true_count<up:
            true_id.append(id)
        else:
            false_id.append(id)
            # #### get the extended true_id

    # ### directly set the transcript_id whose true_count is zeros into the true_id group
    extend_true = list(set(id_qb).difference(set(truth_id)))
    true_id.extend(extend_true)
    # ### concatenate the true and false id
    all_id = list(true_id)
    all_id.extend(false_id)

    # ### add label for the list
    #     set label for every transcript_id(success(true_id,set as 1),fail(false_id,set as 0))
    #     And them we will merge this labeled list with list of properties in order to get a list which include both properties and label of every transcript.
    # add label for the id
    label = []
    for i in range(len(true_id)):
        label.append(1)
    for i in range(len(false_id)):
        label.append(0)

    labeled_id = [all_id,label]
    labeled = list(map(list,zip(*labeled_id)))

    df_labeled_id = pd.DataFrame.from_records(labeled, columns=['Name','label'])
    df_labeled_id.Name = df_labeled_id.Name.astype(str)
    df_labeled_id = df_labeled_id.set_index('Name')
    return df_labeled_id

def cal_penaty(df_labeled_id,seq_id,list_class):
    penaty = 0
    i=0
    for class_member in list_class:
        i = i + 1
        num = class_member[0]
        count = class_member[-1]
        true_count = 0
        for id in class_member[1:-1]:
            true_count += df_labeled_id.loc[seq_id[id]].label
        penaty += (true_count*(num-true_count))*count/num
    return penaty

def loop_for_penaty(root_path):
    df_poly_truth,truth_id = read_poly_truth(root_path)
    df_quant_boot,id_qb = read_bootstraps(root_path)
    # what is use_id : all the id in df_qb_sorted?
    file = open("penaty_2.txt","w")
    df_qb_sorted,sort_qb,use_id = sort_bootstrap(df_quant_boot,truth_id)
    list_class,seq_id = read_eq_classes(root_path)
    # alphaset = np.linspace(0.3,0.475,15)
    alphaset = np.linspace(0.475,0.490,25)
    print(alphaset)
    for alpha in alphaset:
        df_labeled_id = fileter_the_false(df_poly_truth,df_qb_sorted,sort_qb,id_qb,use_id,truth_id,alpha)
        penaty = cal_penaty(df_labeled_id,seq_id,list_class)
        print("alpha = ",alpha)
        print("penaty = ",penaty)
        print("\n")
        file.write(str(alpha)+" "+str(penaty)+"\n")

    file.close()
    return 0

def find_the_min():
    file = open("penaty.txt","r")
    contents = file.readlines()
    min = 100000000
    for line in contents:
        con = line.split(' ')
        penaty = con[1][:-1]
        penaty = float(penaty)
        if penaty<min:
            min = penaty
            min_alpha = float(con[0])
        print(penaty)
    file.close()
    print(min,min_alpha)
    return min,min_alpha
# ## root path
root_path = "../data/poly_mo/"
loop_for_penaty(root_path)
# find_the_min()
