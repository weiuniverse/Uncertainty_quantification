<!-- $theme: gaia -->
### Proj. 3
# Uncertainty Quantification
##### Team members: 
##### Zhengwei Wei, Chong Ou, 
##### Yunqing Yang, Haozhi Qu

---

## Project Description
Salmon is a tool for measuring gene expression, and gives estimates of transcript abundances. Using 'bootstrapping' technique, we can measure the confidence in these estimates. However, this approach tends to underestimate the uncertainty with more transcripts falling out of the interval than expected. 

Our goal is to <b>filter out these failed transcripts, find common properties between them</b>, and come up with a quality score based on those properties which measures our confidence in the estimates. 

---

## Our Approaches
1. Parse the data files (poly_truth.tsv, quant_bootstrap.tsv, quant.sf);
2. Find range of confidence interval;
3. Filter out the failed transcripts;
4. Group data by true/failed transcripts, analyze common properties.
---

## Implementation Details

---

## Parse quant_bootstraps.tsv
This file gives us the bootstrap data (200 rounds of sample taking).  
Note: we are using `DataFrame` from the `pandas` package to easier do data analysis.
```python
quant_bootstraps = tsv.TsvReader(open(root_path + 
				  "quant_bootstraps.tsv"))
quant_boot = [line for line in quant_bootstraps]
```
```python
df_quant_boot = pd.DataFrame.from_records(quant_boot[1:], 
				    columns=quant_boot[0])
```
---

## Parse poly_truth.tsv
This file gives us the true count of each transcript.
```python
# Poly_truth.tsv: true counts for each transcript
poly_truth = open(root_path + "poly_truth.tsv")
lines = poly_truth.readlines()
poly_truth.close()

poly_truth = [['transcript_id','count']]
poly_truth.extend(line[:-1].split('\t') for line in lines)
df_poly_truth = pd.DataFrame.from_records(poly_truth[1:], 
				    columns=poly_truth[0])
```
```python
df_poly_truth['transcript_id']=
		df_poly_truth['transcript_id'].astype(str)
df_poly_truth['count']=df_poly_truth['count'].astype(int)
```

---

## Parse quant.sf
This file gives us some attributes of the transcripts.


```python
quant_file = open(root_path + "quant.sf")
lines = quant_file.readlines()
quant_file.close()
quant = [line[:-1].split('\t') for line in lines]
```
```python
df_quant = pd.DataFrame.from_records(quant[1:], 
					 columns=quant[0])
df_quant.Name = df_quant.Name.astype(str)
df_quant.Length = df_quant.Length.astype(int)
df_quant.EffectiveLength = df_quant.EffectiveLength.astype
						   (float)
df_quant.TPM = df_quant.TPM.astype(float)
df_quant.NumReads = df_quant.NumReads.astype(float)
```

---

We find and retrieve the intersecting transcript ids of poly_truth and quant_bootstraps, and sort each id's data by ascending order. There are transcripts in quant_bootstraps that don't show up in poly_truth, we'll deal with them later.
```python
set_qb_id = set(df_quant_boot.columns)
set_pt_id = set(df_poly_truth.transcript_id)
intersect_ids = set_qb_id & set_pt_id 

sort_qb = []
use_id = []
for id in intersect_ids:
    listed = list(df_quant_boot[id])
    listed.sort()
    use_id.append(id)
    sort_qb.append(listed)
sort_qb = list(map(list,zip(*sort_qb)))
```

---

## Find confidence interval
Since we have already sorted each transcript id's data, we can find an empirical confidence interval of 95% by locating the numbers at index `(total_length) * 2.5%` and `(total_length) * 97.5%`, which would be the lower and upper bound.
```python
df_poly_truth = df_poly_truth.set_index(['transcript_id'])

sum = len(sort_qb)
percent2dot5 = df_qb_sorted.loc[int(sum*0.025)-1]
percent97dot5 = df_qb_sorted.loc[int(sum*0.975)-1]
```

---

## Find the failed transcripts  
Compare the counts given by poly_truth with the lower and upper bound we found earlier. If not in range we treat it as a failed transcript, else true.
```python
true_id = []
false_id = []
for id in use_id:
    down = float(percent2dot5[id])
    up = float(percent97dot5[id])
    true_count = float(df_poly_truth.loc[id])
    if down < true_count < up: 
        true_id.append(id)
    else: 
        false_id.append(id)
```
---

We go back to deal with the 'diff' transcript ids we ignored earlier. The counts of these diff transcript ids are zero, and we assume that these are true transcripts.
```python
true_id.extend(list(set_qb_id.difference(set_pt_id)))
all_id = true_id[:]
all_id.extend(false_id)
```
We add a label of 1 representing true transcripts and 0 representing failed transcripts for easy grouping later on.
```python
label = [1 if i < len(true_id) else 0 for i in 
		      range(len(true_id) + len(false_id))]
labeled_id = [all_id,label]
labeled = list(map(list,zip(*labeled_id)))
```
---

## Common Properties of Failed Transcripts

We group the data by true and failed transcripts, and first observe the mean, std, max and min.  

---

Observing the mean, the average TPM and NumReads of failed transcripts is a lot bigger than the true ones.  

![mean](mean.png)

---

With the std, we find that failed transcripts tend to have a significantly larger variance of TPM and NumReads.  

![std](std.png)

---

Min values don't seem to have much to offer other than maybe a slightly bigger length.

![min](min.png)

---

With max values, the TPM and NumReads attributes seem interesting, which given the info from std, tells us that failed transcripts have a significantly wider range and variance of TPM and NumReads.  

![max](max.png)

---

This looks more clear when we plot the distribution:  

![dis](plot.png)  

---

## Classification Models
We used `sklearn` package for this.  

We use models to predict whether transcripts are true or false. We use data in `poly_mo` folder for our models. After shuffling our data, we make the the properties Length, EffectiveLength, TPM, NumReads as X and label as Y. We then split our data so that we can use 90% of our data as train data and other 10% as test data.  

Note: Because we shuffle our data and randomly select train and test data, the result of our models may be a little bit different whenever you rerun the code.  
```python
sfdata = shuffle(data) # make the data random
input_data = sfdata[['Length','EffectiveLength','TPM','NumReads']]
input_label = sfdata['label']  

train = int(len(input_data)*9/10)
test = train+1
train_data = input_data[0:train]
train_label = input_label[0:train]
test_data = input_data[test:]
test_label = input_label[test:]

def trans_onehot(labels):
    onehot_labels = []
    for label in labels:
        if label==0:
            onehot_labels.append([0,1])
        else:
            onehot_labels.append([1,0])
    return onehot_labels
    
train_onehot_label = trans_onehot(train_label)
test_onehot_label = trans_onehot(test_label)
```

---
**Linear Regression**
First we gave **linear regression** a shot, and we ended up with:  

mean squared error = 0.1550110959836294   
accuracy= 0.6547798066595059  

The results were not ideal. Looks like it  doesn't seem to be  linear, so we tried another classification model.
```
lr = LinearRegression(fit_intercept=True, normalize=False, copy_X=True, n_jobs=1)
lr.fit(dropped_train_data,dropped_train_label)
lr_regr = lr.predict(dropped_test_data)
lr_alpha = 0.5
lr_pred = cal_pred(lr_regr,alpha)

lr_mse = cal_mse(lr_regr,dropped_test_label)
lr_accuracy = cal_accuracy(lr_pred,dropped_test_label)
```

---

**Support Vector Regression**  

With **support vector regression**, we achieved:

mean squared error = 0.06162988735258758   
accuracy = 0.9667338709677419  

96.7% accuracy better than linear regression.  
```python
clf = SVR(C=1.0, epsilon=0.2)
clf.fit(dropped_train_data,dropped_train_label)

svr_regr = clf.predict(dropped_test_data)
svr_mse = cal_mse(svr_regr,dropped_test_label)
svr_alpha = 0.5
svr_pred = cal_pred(svr_regr,svr_alpha)
svr_accuracy = cal_accuracy(svr_pred,dropped_test_label)
```

---  

**Neural network**  

We implemented a simple Nerual network model with 4 hidden layers using Keras. However, this is only a initial model. We will modify the parameters to achieve a comparable accuracy.  

Nerual Network accuracy = 0.7560483870967742  
```python
train_X = dropped_train_data.as_matrix()
train_Y = dropped_train_label.as_matrix()
train_X = np.reshape(train_X,(len(train_X),4))  
train_Y = np.reshape(train_Y,(len(train_Y),1))  
n_samples = train_X.shape[0]

model = Sequential()
model.add(Dense(8, input_shape=(4,)))
model.add(Activation('tanh'))
model.add(Dense(16))
model.add(Activation('tanh'))
model.add(Dense(4))
model.add(Dense(1))

model.compile(optimizer='rmsprop', loss='mse')
model.fit(train_X, train_Y, epochs=100, batch_size=64)

Test = dropped_test_data.as_matrix()
Test_label = dropped_test_label.as_matrix()
nn_regr = model.predict(Test)
nn_pred = cal_pred(nn_regr,0.5)
nn_pred = cal_pred(nn_regr,0.5)
diff = nn_pred-Test_label
diff = diff.tolist()
right_count = diff.count(0)
nn_accuracy = right_count/len(diff)
```

---
## Next Steps
We plan to:  
- Put `eq_classes.txt` into use and adjust the interval  
- Modify parameters in our model to improve accuracy
