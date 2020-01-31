# AveragingCausalHD
R package to estimate causal effects in high dimensions by averaging over many estimators. This implements the approach described in "Averaging causal estimators in high dimensions" by Joseph Antonelli and Matthew Cefalu. The manuscript can be found at

https://arxiv.org/pdf/1906.09303.pdf

To download the R package use the following in R:

```
library(devtools)
install_github(repo = "jantonelli111/AveragingCausalHD")
library(AveragingCausalHD)
```

### Installing other packages
Some of the individual estimators used in the averaging require their own set of R packages. The AveragingCausalHD function will still work if these functions are not installed or loaded, as it will simply drop any estimators that do not have the required packages loaded. Most estimators simply rely on existing CRAN packages such as glmnet, which should install automatically when AveragingCausalHD is loaded. Three estimators, however, rely on their own R packages which are available on Github. If you want to use these three estimators, then you can use the following line of code to install and load them.

```
library(devtools)
install_github(repo = "jantonelli111/DoublyRobustHD")
library(DoublyRobustHD)

install_github(repo = "jantonelli111/HDconfounding")
library(HDconfounding)

install_github("swager/balanceHD")
library(balanceHD)
```

### How to use the software

The software estimates the average treatment effect of a binary treatment on a continuous outcome while adjusting for a potentially high-dimensional set of covariates. The software has 10 built in estimators for estimating the treatment effect, and then combines the individual estimators to provide a more robust estimate of the treatment effect. First, we need to simulate data for this scenario

```{r, eval=FALSE}
n <- 100
p <- 100

beta.c <- c(0.75,1, 0.6, -0.8, -0.7,rep(0, p-5))

gamma <- c(0.15,0.2,0, 0, -0.4, rep(0, p-5))

beta <- 1

sigma <- matrix(0.3, p,p)
diag(sigma) <- 1

x <- mvtnorm::rmvnorm(n, sigma=sigma)
t <- as.numeric((x %*% gamma + rnorm(n)) > 0)
y <- 0 + t + x %*% beta.c + rnorm(n, sd=1)
```

And now that we have the data, we can apply the main function as follows:

```{r, eval=FALSE}
results <- AveragingCausalHD(y=y, t=t, x=x)
```

The main estimator results can be found with the following command
```{r, eval=FALSE}
results$averaged
```

If the individual estimator results are also of interest, then the following command will show them

```{r, eval=FALSE}
results$individual_estimators
```

### Specifying which estimators to use

If not specified, then the averaged estimator will use all estimators that do not rely on MCMC to save computation time. If the full set of 10 estimators is desired, this can be specified as

```{r, eval=FALSE}
AveragingCausalHD(y=y, t=t, x=x,
                estimators = c("DoublePS", "DRlasso", "Debiasing",
                               "DML", "DMLpost_selection", "TMLElasso",
                               "TMLEscreen","HDmatching", "HDbayes", "HDC"))
```

The HDbayes and HDC estimators might take longer to fit as they are Bayesian approaches, so they can be removed if computation time is a concern. 

### Adding user-specified estimators

If the user has their own estimates of the treatment effect that they would like to enter into the averaging estimator, this can be done using the AdditionalEstimates and AdditionalSEs arguments. Both must be specified in order for this to work. Let's suppose that I had an estimate of the treatment effect of 1.1 with a corresponding standard error estimate of 0.25. I could couple this with the estimators already in AveragingCausalHD as follows:


```{r, eval=FALSE}
AveragingCausalHD(y=y, t=t, x=x,
                  AdditionalEstimates = 1.1,
                  AdditionalSEs = 0.25)
```
