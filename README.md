# R package YoudenDRM
This package is developed based on [Yuan et al. (2021)](https://onlinelibrary.wiley.com/doi/abs/10.1002/cjs.11600). 
This package provides the estimators and confidence intervals for the Youden index and optimal cutoff point under the density ratio models.


### Table of Contents
**[Installation](#installation)**<br>
**[Functions](#functions)**<br>
**[Usage](#usage)**<br>
**[References](#references)**<br>
## Installation
Install this package from Github with 
```r
# If the package "devtool" has not been installed, install it fisrt with 
#install.packages("devtools")
library(devtools)
devtools::install_github("MengYuan-UW/YoudenDRM")
library("YoudenDRM")
```
## Functions
This package contains the following functions:
- The Density Ratio Model (`DRM`): a function to fit the DRM.
- The inference on the Youden index and cutoff point (`Youden`): a function to estimate the Youden index and cutoff point as well as construct their confidence intervals.
- The goodnees-of-fit test for the DRM (`goodnessFit`): a test to check the validity of the DRM with a pre-specified basis function.


## Usage
We provide two examples.
- Example 1: data without lower limit of detection
```r
library("YoudenDRM")
# generate the fully observed data
set.seed(123456)
x = rlnorm(50, meanlog = 2.5, sdlog = sqrt(0.09))
y = rlnorm(50, meanlog = 2.87,sdlog = sqrt(0.25))
# basis function
qt = c("log(t)","log(t)^2")
# perform the goodness-of-fit test for the density ratio model with basis function
set.seed(123456)
goodnessFit(x,y,qt,B = 1000)
# obtain the estimate, asymptotic standard deviation (ASD), and the 95% confidence intervals (lower bound and uppper bound) of the Youden index and optimal cutoff point
Youden(x,y,qt,CItype ="logit-DRM")
```
- Example 2: data with a lower limit of detection
```r
# create the LLOD
r = qlnorm(0.15,meanlog = 2.5, sdlog = sqrt(0.09))
# get the data with the lower limit of detection r
x = x[x>=r]
y = y[y>=r]
# obtain the estimate, asymptotic standard deviation (ASD), and the 95% confidence intervals (lower bound and uppper bound) of the Youden index and optimal cutoff point
Youden(x,y,qt,r,totalSize = c(50,50),CItype ="logit-DRM")
```

## References
Yuan M, Li P, Wu C (2021). “Semiparametric inference of the Youden index and the optimal cut-off point under density ratio models.” The Canadian Journal of Statistics, 49, 965 - 986. [PDF](https://onlinelibrary.wiley.com/doi/abs/10.1002/cjs.11600)

