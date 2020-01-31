#' Estimate the average treatment effect while adjusting for a high dimensional set of covariates
#'
#'
#' @param y                            The outcome to be analyzed
#' @param t                            The treatment to be analyzed
#' @param x                            An n by p matrix of covariates to adjust for
#' @param estimators                   A list of estimators to be included in the averaging estimator.
#'                                     The full list of estimators is given by c("DoublePS", "DRlasso", "Debiasing",
#'                                     "DML", "DMLpost_selection", "TMLElasso", "TMLEscreen", "HDC", "HDbayes", "HDmatching")
#'                                     The HDbayes and HDC estimators are Bayesian and may take time to
#'                                     run on large data sets, however, the rest should run quickly. It is
#'                                     recommended to use as many of these estimators in the averaging as possible
#'
#' @param AdditionalEstimates          A vector of additional point estimates the user can include in the averaging.
#'                                     This should only be used if the user has fit other estimators that are not
#'                                     already available in the list above, and would like to include them in the 
#'                                     averaging.
#'                                     
#' @param AdditionalSEs                A vector of additional standard errors that correspond to the point estimates
#'                                     from the estimators making up AdditionalEstimates. 
#'                                    
#' @param trim                         An indicator of whether the averaging should remove the highest and lowest 
#'                                     point estimates from the averaging procedure          
#'                                                                                                                                                                             
#' @return A list with the estimates and standard errors from each individual estimator, as well as
#'         from the averaged estimator, which combines all of the individual estimators
#'
#' @export
#' @examples
#'
#' n <- 100
#' p <- 100
#' 
#' beta.c <- c(0.75,1, 0.6, -0.8, -0.7,rep(0, p-5))
#' 
#' gamma <- c(0.15,0.2,0, 0, -0.4, rep(0, p-5))
#' 
#' beta <- 1
#' 
#' sigma <- matrix(0.3, p,p)
#' diag(sigma) <- 1
#' 
#' x <- mvtnorm::rmvnorm(n, sigma=sigma)
#' t <- as.numeric((x %*% gamma + rnorm(n)) > 0)
#' y <- 0 + t + x %*% beta.c + rnorm(n, sd=1)
#' 
#' AveragingCausalHD(y=y, t=t, x=x)

AveragingCausalHD = function(y,t,x,estimators = c("DoublePS", "DRlasso", "Debiasing",
                                    "DML", "DMLpost_selection", "TMLElasso",
                                    "TMLEscreen", "HDmatching"),
                     AdditionalEstimates = NULL, AdditionalSEs = NULL,
                     trim = FALSE) {
 
  ## Make sure that all required packages for estimators are loaded
  packageList = (.packages())
  
  ## remove estimators that don't have correct packages
  if ("DoublePS" %in% estimators) {
    if (!"glmnet" %in% packageList) {
      estimators = estimators[-which(estimators == "DoublePS")]
      print("glmnet needs to be loaded for the Double PS estimator")
      print("Therefore the Double PS estimator has been dropped")
    }
  }
  
  if ("DRlasso" %in% estimators) {
    if (!"glmnet" %in% packageList) {
      estimators = estimators[-which(estimators == "DRlasso")]
      print("glmnet needs to be loaded for the DRlasso estimator")
      print("Therefore the DRlasso estimator has been dropped")
    } else {
      library(glmnet)
    }
  }
  
  if ("Debiasing" %in% estimators) {
    if (!"glmnet" %in% packageList | ! "balanceHD" %in% packageList) {
      estimators = estimators[-which(estimators == "Debiasing")]
      print("glmnet and balanceHD need to be loaded for the Debiasing estimator")
      print("Therefore the Debiasing estimator has been dropped")
    } else {
      library(glmnet)
      library(balanceHD)
    }
  }
  
  if ("DML" %in% estimators) {
    if (!"glmnet" %in% packageList) {
      estimators = estimators[-which(estimators == "DML")]
      print("glmnet needs to be loaded for the DML estimator")
      print("Therefore the DML estimator has been dropped")
    } else {
      library(glmnet)
    }
  }
  
  if ("DMLpost_selection" %in% estimators) {
    if (!"glmnet" %in% packageList) {
      estimators = estimators[-which(estimators == "DMLpost_selection")]
      print("glmnet needs to be loaded for the DMLpost_selection estimator")
      print("Therefore the DMLpost_selection estimator has been dropped")
    } else {
      library(glmnet)
    }
  }
  
  if ("TMLElasso" %in% estimators) {
    if (!"glmnet" %in% packageList | !"tmle" %in% packageList) {
      estimators = estimators[-which(estimators == "TMLElasso")]
      print("glmnet and tmle need to be loaded for the TMLElasso estimator")
      print("Therefore the TMLElasso estimator has been dropped")
    } else {
      library(glmnet)
      library(tmle)
    }
  }
  
  if ("TMLEscreen" %in% estimators) {
    if (!"glmnet" %in% packageList | !"tmle" %in% packageList) {
      estimators = estimators[-which(estimators == "TMLEscreen")]
      print("glmnet and tmle need to be loaded for the TMLEscreen estimator")
      print("Therefore the TMLEscreen estimator has been dropped")
    } else {
      library(glmnet)
      library(tmle)
    }
  }
  
  if ("HDC" %in% estimators) {
    if (!"HDconfounding" %in% packageList) {
      estimators = estimators[-which(estimators == "HDC")]
      print("HDconfounding needs to be loaded for the HDC estimator")
      print("Therefore the HDC estimator has been dropped")
    } else {
      library(HDconfounding)
    }
  }
  
  if ("HDbayes" %in% estimators) {
    if (!"DoublyRobustHD" %in% packageList) {
      estimators = estimators[-which(estimators == "HDbayes")]
      print("DoublyRobustHD needs to be loaded for the HDbayes estimator")
      print("Therefore the HDbayes estimator has been dropped")
    } else {
      library(DoublyRobustHD)
    }
  }
  
  if ("HDmatching" %in% estimators) {
    if (!"Matching" %in% packageList | !"glmnet" %in% packageList) {
      estimators = estimators[-which(estimators == "HDmatching")]
      print("Matching and glmnet need to be loaded for the HDmatching estimator")
      print("Therefore the HDmatching estimator has been dropped")
    } else {
      library(Matching)
      library(glmnet)
    }
  }
  
  results = data.frame(estimator = estimators,
                       est = rep(NA, length(estimators)),
                       se = rep(NA, length(estimators)),
                       CIlower = rep(NA, length(estimators)),
                       CIupper = rep(NA, length(estimators)))
  
  counter = 1
  for (estimator in estimators) {
    if (estimator == "DoublePS") {
      print("Fitting DoublePS estimator now")
      tempResults = NULL
      tempResults = DoublePS(y=y, t=t, x=x)
    } else if (estimator == "DRlasso") {
      print("Fitting DRlasso estimator now")
      tempResults = NULL
      tempResults = DRlasso(y=y, t=t, x=x)
    } else if (estimator == "Debiasing") {
      print("Fitting Debiasing estimator now")
      tempResults = NULL
      tempResults = Debiasing(y=y, t=t, x=x)
    } else if (estimator == "DML") {
      print("Fitting DML estimator now")
      tempResults = NULL
      tempResults = DML(y=y, t=t, x=x)
    } else if (estimator == "DMLpost_selection") {
      print("Fitting DMLpost_selection estimator now")
      tempResults = NULL
      tempResults = DMLps(y=y, t=t, x=x)
    } else if (estimator == "TMLElasso") {
      print("Fitting TMLElasso estimator now")
      tempResults = NULL
      tempResults = TMLElasso(y=y, t=t, x=x)
    } else if (estimator == "TMLEscreen") {
      print("Fitting TMLEscreen estimator now")
      tempResults = NULL
      tempResults = TMLEscreen(y=y, t=t, x=x)
    } else if (estimator == "HDC") {
      print("Fitting HDC estimator now")
      tempResults = NULL
      tempResults = HDC(y=y, t=t, x=x)
    } else if (estimator == "HDbayes") {
      print("Fitting HDbayes estimator now")
      tempResults = NULL
      tempResults = DRB(y=y, t=t, x=x)
    } else if (estimator == "HDmatching") {
      print("Fitting HDmatching estimator now")
      tempResults = NULL
      tempResults = DRmatch(y=y, t=t, x=x)
    }
    
    results$est[counter] = tempResults$est
    results$se[counter] = (tempResults$Interval[2] - tempResults$Interval[1])/3.92
    results$CIlower[counter] = tempResults$Interval[1]
    results$CIupper[counter] = tempResults$Interval[2]
    
    counter = counter + 1
  }

  
  ## add in results from user specified estimators if they exist
  if (is.null(AdditionalEstimates) == FALSE &
      is.null(AdditionalSEs) == FALSE) {
    addResults = data.frame(estimator = paste("User estimate", 1 : length(AdditionalSEs)),
                            est = AdditionalEstimates,
                            se = AdditionalSEs,
                            CIlower = AdditionalEstimates - 1.96*AdditionalSEs,
                            CIupper = AdditionalEstimates + 1.96*AdditionalSEs)

    results = rbind(results, addResults)
  }
  
  individualResults = results
  
  ## First trim results if that option is set to true
  if (trim == TRUE & nrow(results) < 4) {
    print("Trimming only applies with 4 or more estimators")
    print("Analysis will proceed without trimming")
    trim = FALSE
  }
  
  if (trim == TRUE) {
    wMax = which.max(results$est)
    wMin = which.min(results$est)
    results = results[-c(wMin, wMax),]
  }
  
  ## Now remove any NAs for the estimation
  wNA = which(!is.na(results$est))
  
  ## error if there are no estimators left
  if (length(wNA) == 0) {
    print("Error: No individual estimators are left to average over")
  }
  
  results = results[wNA,]
  
  ## Now perform averaging on the remaining estimators
  varEst = 0
  for (i in 1 : nrow(results)) {
    for (j in 1 : nrow(results)) {
      varEst = varEst + results$se[i]*results$se[j]
    }
  }
  
  varEst = varEst / (nrow(results)^2)
  
  avgEst = data.frame(est = mean(results$est),
                      se = sqrt(varEst),
                      CIlower = mean(results$est) - 1.96*sqrt(varEst),
                      CIupper = mean(results$est) + 1.96*sqrt(varEst))
  
  l = list(averaged = avgEst,
           individual_estimators = individualResults)
  
  return(l)
}
  
