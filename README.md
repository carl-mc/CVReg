# CVReg: A Covoting Regression Model and Associated Methods

## Introduction

This package helps researchers to prepare data and estimate the covoting
regression model (CVR) developed by Müller-Crepon and Bormann (APSR,
forthcoming) and briefly outlined below.  
When using the method and/or package, please cite:

Müller-Crepon, C., & Bormann, N.-C. (2025). Coethnics Covote in Africa:
Studying Electoral Cleavages with a Covoting Regression Model. *American
Political Science Review*, forthcoming.

You can find the publication for all technical details [here]() and as
ungated version
[here](http://www.carlmueller-crepon.org/publication/ethnic_parties/).

## The method

## Installation

You can directly download and install the CVReg package from GitHub.
Upon installation, the package should automatically install all
necessary R dependencies.

``` r
library(devtools)
# install_github(repo = "carl-mc/CVReg")
devtools::load_all(".")
```

## Preparing your data

Users of the CVR model get started by transforming their data (e.g., a
survey) into pairwise observations.

``` r
library(CVReg)
# Load some synthetic data
data("synthvote")

# Generate pairs
pair.dt <- data_to_pairs(df = synthvote, structure = "full")

# Inspect
head(pair.dt)
```

    ##    a.resp a.group  a.income a.unit       a.x        a.y   a.prtvt  a.ID b.resp
    ##     <int>  <char>     <num>  <num>     <num>      <num>    <char> <int>  <int>
    ## 1:      1       A 0.5541771      1 0.5308088 0.07287736 MarsParty     1      2
    ## 2:      1       A 0.5541771      1 0.5308088 0.07287736 MarsParty     1      3
    ## 3:      2       A 0.6882752      1 0.6848609 0.87430875 MarsParty     2      3
    ## 4:      1       A 0.5541771      1 0.5308088 0.07287736 MarsParty     1      4
    ## 5:      2       A 0.6882752      1 0.6848609 0.87430875 MarsParty     2      4
    ## 6:      3       A 0.6580576      1 0.3832834 0.19887281 MarsParty     3      4
    ##    b.group  b.income b.unit       b.x       b.y   b.prtvt  b.ID
    ##     <char>     <num>  <num>     <num>     <num>    <char> <int>
    ## 1:       A 0.6882752      1 0.6848609 0.8743088 MarsParty     2
    ## 2:       A 0.6580576      1 0.3832834 0.1988728 MarsParty     3
    ## 3:       A 0.6580576      1 0.3832834 0.1988728 MarsParty     3
    ## 4:       B 0.6633427      1 0.9549880 0.5345898 MoonParty     4
    ## 5:       B 0.6633427      1 0.9549880 0.5345898 MoonParty     4
    ## 6:       B 0.6633427      1 0.9549880 0.5345898 MoonParty     4

Based on the pairs of observations a and b, we can now create a series
of variables that record similarities and shared characteristics in each
pair.

``` r
# Process pair level data to create pairwise similarity/difference measures
pair.dt <- process_pairs(df = pair.dt,
                         vars =  c("prtvt.eq", "group.eq", "income.sim"))

# Inspect
summary(subset(pair.dt,,c("prtvt.eq", "group.eq", "income.sim")))
```

    ##     prtvt.eq       group.eq        income.sim        
    ##  Min.   :0.00   Min.   :0.0000   Min.   :-0.9986159  
    ##  1st Qu.:0.00   1st Qu.:0.0000   1st Qu.:-0.5085070  
    ##  Median :0.00   Median :1.0000   Median :-0.2982530  
    ##  Mean   :0.34   Mean   :0.5298   Mean   :-0.3392353  
    ##  3rd Qu.:1.00   3rd Qu.:1.0000   3rd Qu.:-0.1370884  
    ##  Max.   :1.00   Max.   :1.0000   Max.   :-0.0000033

## Fitting a CVR

Using the pair level data on covoting, shared group membership and
income similarity, we can now run a simple linear regression model.

``` r
library(fixest)
lm <- feols(prtvt.eq ~ group.eq + income.sim, data = pair.dt,
            cluster = c("a.ID", "a.ID"))
summary(lm)
```

    ## OLS estimation, Dep. Var.: prtvt.eq
    ## Observations: 124,750
    ## Standard-errors: Clustered (a.ID) 
    ##             Estimate Std. Error  t value   Pr(>|t|)    
    ## (Intercept) 0.294947   0.009317 31.65674  < 2.2e-16 ***
    ## group.eq    0.146161   0.015072  9.69722  < 2.2e-16 ***
    ## income.sim  0.095371   0.017190  5.54809 4.6979e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## RMSE: 0.46752   Adj. R2: 0.025988

Or we can use `the estimate_cvreg()` function, which implements some
additional methods, such as computing conditional within- and
between-group covoting probabilities, or clustered standard errors for
dyadic data.

``` r
cvr <- estimate_cvreg(prtvt.eq ~ group.eq + income.sim, data = pair.dt,
                      cvr.det.var = "group.eq", cluster = c("a.ID", "a.ID"))
all(cvr$coefficient == lm$coefficient)
```

    ## [1] TRUE

``` r
summary(cvr)
```

    ## OLS estimation, Dep. Var.: prtvt.eq
    ## Observations: 124,750
    ## Weights: weight.vec
    ## Standard-errors: Clustered (a.ID) 
    ##             Estimate Std. Error  t value   Pr(>|t|)    
    ## (Intercept) 0.294947   0.009317 31.65674  < 2.2e-16 ***
    ## group.eq    0.146161   0.015072  9.69722  < 2.2e-16 ***
    ## income.sim  0.095371   0.017190  5.54809 4.6979e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## RMSE: 0.46752   Adj. R2: 0.025988

``` r
unlist(cvr$sum.ls[c(1,2)])
```

    ## covote_rate  grpeq_mean 
    ##   0.3400321   0.5298116

``` r
do.call(cbind, cvr$sum.ls[-c(1,2)])
```

    ##     covoting_within_grp covoting_between_grp covoting_ratio sharedgrp_covoters
    ## est         0.408755310          0.262594220       1.556604         0.63689086
    ## se          0.009404921          0.006857815       0.070748         0.01465402
    ##     sharedgrp_noncovoters sharedgrp_ratio
    ## est           0.474641709      1.34183501
    ## se            0.007550119      0.01899355

## Contrast with conditional covoting probabilities from multinomial model

``` r
est_covote_indlevel(out = "prtvt", cov = c("group", "income"), data = synthvote, method = "ols")
```

    ##                 beta  covoting_within_grp covoting_between_grp 
    ##            0.1426601            0.4083003            0.2656402

## Add cross-cleavage interaction terms for improved inference

``` r
m.int <- estimate_cvreg(prtvt.eq ~ group.eq + income.sim, data = pair.dt, 
                      cluster = c("a.ID", "a.ID"), main.int.var = "group.eq",
                      int.vars = c("income.sim"))
```

    ## [1] "LASSO Interaction Selection"

``` r
summary(m.int)
```

    ## OLS estimation, Dep. Var.: prtvt.eq
    ## Observations: 124,750
    ## Weights: weight.vec
    ## Standard-errors: Clustered (a.ID) 
    ##                Estimate Std. Error   t value   Pr(>|t|)    
    ## (Intercept)    0.264375   0.008472 31.203931  < 2.2e-16 ***
    ## group.eq       0.140631   0.014721  9.552916  < 2.2e-16 ***
    ## income.sim    -0.003280   0.015774 -0.207957 8.3535e-01    
    ## int__c1_c1     0.149837   0.025981  5.767088 1.4172e-08 ***
    ## int__group2.2 -0.025566   0.015574 -1.641586 1.0131e-01    
    ## int__m1_c1    -0.096686   0.016314 -5.926747 5.7752e-09 ***
    ## int__m2_c1     0.066747   0.013221  5.048661 6.2509e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## RMSE: 0.465822   Adj. R2: 0.033017

## Feedback, comments, questions

We are very grateful for any bug reports, feedback, questions, or
contributions to this package. Please report any issues here or write to
c.a.muller-crepon \[at\] lse.ac.uk .
