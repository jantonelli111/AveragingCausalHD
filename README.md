# AveragingCausalHD
R package to estimate causal effects in high dimensions by averaging over many estimators. This implements the approach described in "Averaging causal estimators in high dimensions" by Joseph Antonelli and Matthew Cefalu. The manuscript can be found at

https://arxiv.org/pdf/1906.09303.pdf

To download the R package use the following in R:

```
library(devtools)
install_github(repo = "jantonelli111/AveragingCausalHD")
library(AveragingCausalHD)
```

# How to use the software

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
results = AveragingCausalHD(y=y, t=t, x=x)
```
