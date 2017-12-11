# Abstract

# Introduction
## Desribe the problem


# Methods
## Goal
The Goal of this project is to modify the bootstraping value close to true count. Therefore,
we should predict the difference between bootstraping mean value and true count. We set $$ distance = mean - true_count$$. And we have some features to predict the distance, which means
distance = f(features). The features here include the properties of the transcript and the mean value and standard deviation value of bootstraping. Therefore our task is to use statistic or machine learning methods to approximate the function distance = f(features).

## Machine Learning Model
We tried several models to predict distance.
### Linear Regression
This model is used to fit linear model. We first tried this model to get a result.
Given a dataset $${distance_i,X_i}_i^n, X is features$$,
The model take the form: $$distance=aX+b $$.
