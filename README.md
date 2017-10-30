# Uncertainty_quantification
## some useful docs:
https://zhuanlan.zhihu.com/p/26319993

## data description：
You are provided with two sample bootstrapping outputs each containing the four following files:

Quant_bootstraps.tsv  containing the matrix of bootstrap experiments containing the final count for each transcript in each round of bootstrapping with a row be a bootstrap output and columns be list of transcripts.

Poly_truth.tsv: true counts for each transcript

Eq_classes.txt: list of equivalence classes and their information

Quant.sf: estimated quantifications for each transcript

## Project Description:

Salmon is a state-of-the-art tool for measuring gene expression (quantifying the abundance of different RNA transcripts in an experiment). One of the primary benefits of this method is that it is orders of magnitude faster than the competition, while producing results of the same or better accuracy. Salmon determines expression levels by solving a maximum likelihood problem. A result of this formulation is that one often gets accurate estimates of the transcript abundances, but has no notion of confidence in these predictions. That is, predictions can be highly specific and highly accurate, or highly specific and highly inaccurate --- this depends on the “shape” of the likelihood function, and how optimization proceeds. Thus, it is very valuable to provide some measure of confidence in the estimates that are inferred. There are a number of ways to estimate such confidence. One way is “bootstrapping”. This approach treats the observed sample data as the population, samples from the original data a number of times, and reruns the maximum likelihood estimator independently on all of these samples. By looking at the distribution of these different runs, one can form an empirical confidence interval, that provides a notion of the uncertainty of the maximum likelihood point estimate. However, a close inspection on simulated data demonstrates that these empirical intervals fail to capture the uncertainty adequately (i.e., they tend to underestimate the uncertainty).  

The goal here is to determine what are the specific properties, if any, of the transcripts that tend to fall outside the predicted uncertainty interval. To be exact, we define the problem as follows: for each transcript expression, if the truth value is not within the 95% confidence interval  (i.e. greater than at least 2.5% of the bootstrap samples and less than the upper 97.5 percentile of the bootstrap samples) we count the transcript as a failure. We expect that, in statistically valid cases, only 5% of transcripts fall outside of the 95% confidence intervals. The steps of the project would be to filter the failing transcripts and then find some common properties between them, using sequence similarity or transcript quantification estimates or other variables output by Salmon during the bootstrapping phase. The second part is to come up with a quality score based on these properties. This score will reflect the level of confidence we have in the transcript-level estimates.

**Inputs:**   

The following information about each transcript:
- Length  
- Effective length
- Set of equivalence classes + their counts and transcript weight distributions
- Count of ambiguous vs uniquely mapped reads
- Transcripts tpm distributions
- True transcript tpm

**Output:**  

- Property values + A Quality Measurement of that property

**Midway Result:**  

Filter the faulty transcripts and provide some statistics regarding their properties.

We will create a project item on blackboard that you should upload the code and one page report (short explanation of your approach and your results) there until Nov. 6th.

We will also create a google sheet and schedule the 10 minute slots on Tuesday Nov. 7th that you are supposed to present your work. I'll send the link to Piazza and you can put your group's names there. Each group will have 7 minutes to explain the code, their progress and show the midway results.
