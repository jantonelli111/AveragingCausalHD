DoublePS = function(y, t, x) {
  ## Double selection method with 1se CV
  
  time1 = Sys.time()
  fit.lasso <- cv.glmnet(cbind(t, x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  fit.lasso.treat <- cv.glmnet(x=as.matrix(x), y=t, intercept=TRUE, family="binomial")
  
  activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
  activeY <- which(coef(fit.lasso, s='lambda.1se')[-c(1,2)] != 0)
  DoubleActive = union(activeX, activeY)
  
  if (length(DoubleActive) == 0) {
    ModelDouble = lm(y ~ t)
  } else {
    ModelDouble = lm(y ~ t + x[,DoubleActive])
  }
  
  est = ModelDouble$coef[2]
  DoubleIntervalAsymp = c(est - 1.96*sqrt(vcov(ModelDouble)[2,2]),
                          est + 1.96*sqrt(vcov(ModelDouble)[2,2]))
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = DoubleIntervalAsymp,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

DRlasso = function(y, t, x) {
  time1 = Sys.time()
  z = t
  fit.x <- cv.glmnet(x, z, intercept = TRUE, family = 'binomial')
  fit.y2 <- cv.glmnet(cbind(z,x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  
  active.x <- which(coef(fit.x)[-1] != 0)
  active.y <- which(coef(fit.y2)[-c(1,2)] != 0)
  
  if (length(active.x) == 0) {
    ps.mod.farrell = glm(z ~ 1, family='binomial')
  } else {
    ps.mod.farrell = glm(z ~ x[,active.x], family='binomial')
  }
  
  if (length(active.y) == 0) {
    out.mod.farrell = lm(y ~ z)
  } else {
    out.mod.farrell = lm(y ~ z + x[,active.y]) 
  }
  
  ps1.farrell = ps.mod.farrell$fitted.values
  ps0.farrell = 1 - ps1.farrell
  
  out1.farrell = cbind(rep(1,n), rep(1,n), x[,active.y]) %*% out.mod.farrell$coefficients
  out0.farrell = cbind(rep(1,n), rep(0,n), x[,active.y]) %*% out.mod.farrell$coefficients
  
  est = sum(z*(y - out1.farrell)/ps1.farrell + out1.farrell)/n -
    sum((1-z)*(y - out0.farrell)/ps0.farrell + out0.farrell)/n
  
  part1 = mean((z*(y - out1.farrell)^2)/(ps1.farrell^2))
  part2 = mean((out1.farrell - mean(y[z==1]))^2)
  var1 = part1 + part2
  
  part2 = mean((out1.farrell - mean(y[z==1])) * (out0.farrell - mean(y[z==0])))
  cov01 = part2
  
  part1 = mean(((1-z)*(y - out0.farrell)^2)/(ps0.farrell^2))
  part2 = mean((out0.farrell - mean(y[z==0]))^2)
  var0 = part1 + part2
  
  CovMat = matrix(c(var0, cov01, cov01, var1), nrow=2)/n
  
  totalVar = c(-1,1) %*% CovMat %*% c(-1,1)
  totalSE = sqrt(totalVar)
  
  time2 = Sys.time()
  
  Interval = c(est - 1.96*totalSE, est + 1.96*totalSE)
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

Debiasing = function(y, t, x) {
  time1 = Sys.time()
  tau.hat = residualBalance.ate(x, y, t, estimate.se = TRUE)
  est = tau.hat[1]
  Interval = c(tau.hat[1] - 1.96*tau.hat[2],
               tau.hat[1] + 1.96*tau.hat[2])
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

DML = function(y, t, x) {
  time1 = Sys.time()
  n = length(y)
  K = 5
  est = rep(NA, K)
  psiTranspose = rep(NA, K)
  psiA = rep(NA, K)
  
  for (k in 1 : 5) {
    obs = ((k-1)*(n/K) + 1) : (k*(n/K))
    
    yc = y[-obs]
    tc = t[-obs]
    xc = x[-obs,]
    
    yi = y[obs]
    ti = t[obs]
    xi = x[obs,]

    set.seed(k)
    
    fitY1 <- cv.glmnet(cbind(tc,xc), yc, intercept = TRUE, 
                       penalty.factor = rep(c(0,1), c(1, p)))
    
    fitT1 <- cv.glmnet(cbind(xc), tc, intercept = TRUE, family="binomial")
    
    g0hat = as.vector(cbind(rep(1,length(obs)), rep(0,length(obs)), xi) %*% coef(fitY1))
    
    d0hat = as.vector(predict(fitT1, newx = cbind(xi), type="response"))

    est[k] = mean((yi - g0hat)*(ti - d0hat)) / mean(ti*(ti - d0hat))
    
    psiA[k] = mean(-ti*(ti - d0hat))
  }
  
  theta0 = mean(est)

  for (k in 1 : 5) {
    obs = ((k-1)*(n/K) + 1) : (k*(n/K))
    
    yc = y[-obs]
    tc = t[-obs]
    xc = x[-obs,]
    
    yi = y[obs]
    ti = t[obs]
    xi = x[obs,]
    
    set.seed(k)
    
    fitY1 <- cv.glmnet(cbind(tc,xc), yc, intercept = TRUE, 
                       penalty.factor = rep(c(0,1), c(1, p)))
    
    fitT1 <- cv.glmnet(cbind(xc), tc, intercept = TRUE, family="binomial")
    
    g0hat = as.vector(cbind(rep(1,length(obs)), rep(0,length(obs)), xi) %*% coef(fitY1))
    
    d0hat = as.vector(predict(fitT1, newx = cbind(xi), type="response"))
    
    psiTranspose[k] = mean(((yi - ti*theta0 - g0hat)*(ti - d0hat))^2)
    
  }
  
  J0 = mean(psiA)
  varEst = (1/J0)^2 * mean(psiTranspose) / n
  sdEst = sqrt(varEst)
  
  time2 = Sys.time()
  
  l = list(est = theta0,
           Interval = c(theta0 - 1.96*sdEst,
                        theta0 + 1.96*sdEst),
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
  
}

DMLps = function(y, t, x) {
  time1 = Sys.time()
  K = 5
  est = rep(NA, K)
  psiTranspose = rep(NA, K)
  psiA = rep(NA, K)
  
  for (k in 1 : 5) {
    obs = ((k-1)*(n/K) + 1) : (k*(n/K))
    
    yc = y[-obs]
    tc = t[-obs]
    xc = x[-obs,]
    
    yi = y[obs]
    ti = t[obs]
    xi = x[obs,]
    
    set.seed(k)
    
    fitY1 <- cv.glmnet(cbind(tc,xc), yc, intercept = TRUE, 
                       penalty.factor = rep(c(0,1), c(1, p)))
    
    fitT1 <- cv.glmnet(cbind(xc), tc, intercept = TRUE, family="binomial")
    
    wy = which(coef(fitY1)[-c(1,2)] != 0)
    wt = which(coef(fitT1)[-c(1)] != 0)
    
    modY1 = lm(yc ~ tc)
    g0hat = cbind(rep(1,length(obs)), rep(0,length(obs))) %*% coef(modY1)
    if (length(wy) > 0) {
      modY1 = lm(yc ~ tc + xc[,wy])
      g0hat = cbind(rep(1,length(obs)), rep(0,length(obs)), xi[,wy]) %*% coef(modY1)
    }
    
    modT1 = glm(tc ~ 1, family=binomial)
    d0hat = expit(rep(1,length(obs)) * coef(modT1))
    if (length(wt) > 0) {
      modT1 = glm(tc ~ xc[,wt], family=binomial)
      d0hat = expit(cbind(rep(1,length(obs)), xi[,wt]) %*% coef(modT1))
    }
    
    est[k] = mean((yi - g0hat)*(ti - d0hat)) / mean(ti*(ti - d0hat))
    
    psiA[k] = mean(-ti*(ti - d0hat))
  }
  
  theta0 = mean(est)
  
  for (k in 1 : 5) {
    obs = ((k-1)*(n/K) + 1) : (k*(n/K))
    
    yc = y[-obs]
    tc = t[-obs]
    xc = x[-obs,]
    
    yi = y[obs]
    ti = t[obs]
    xi = x[obs,]
    
    set.seed(k)
    
    fitY1 <- cv.glmnet(cbind(tc,xc), yc, intercept = TRUE, 
                       penalty.factor = rep(c(0,1), c(1, p)))
    
    fitT1 <- cv.glmnet(cbind(xc), tc, intercept = TRUE, family="binomial")
    
    wy = which(coef(fitY1)[-c(1,2)] != 0)
    wt = which(coef(fitT1)[-c(1)] != 0)
    
    modY1 = lm(yc ~ tc)
    g0hat = cbind(rep(1,length(obs)), rep(0,length(obs))) %*% coef(modY1)
    if (length(wy) > 0) {
      modY1 = lm(yc ~ tc + xc[,wy])
      g0hat = cbind(rep(1,length(obs)), rep(0,length(obs)), xi[,wy]) %*% coef(modY1)
    }
    
    modT1 = glm(tc ~ 1, family=binomial)
    d0hat = expit(rep(1,length(obs)) * coef(modT1))
    if (length(wt) > 0) {
      modT1 = glm(tc ~ xc[,wt], family=binomial)
      d0hat = expit(cbind(rep(1,length(obs)), xi[,wt]) %*% coef(modT1))
    }
    
    psiTranspose[k] = mean(((yi - ti*theta0 - g0hat)*(ti - d0hat))^2)
    
  }
  
  J0 = mean(psiA)
  varEst = (1/J0)^2 * mean(psiTranspose) / n
  sdEst = sqrt(varEst)
  
  time2 = Sys.time()
  
  l = list(est = theta0,
           Interval = c(theta0 - 1.96*sdEst,
                        theta0 + 1.96*sdEst),
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
  
}

TMLElasso = function(y, t, x) {
  time1 = Sys.time()
  result1 <- tmle(Y=y,A=t,W=x,
                  g.SL.library = c("SL.glmnet"),
                  Q.SL.library = c("SL.glmnet"))
  
  est = result1$estimates$ATE$psi
  Interval = c(result1$estimates$ATE$CI[1],
               result1$estimates$ATE$CI[2])
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

TMLEscreen = function(y, t, x) {
  
  time1 = Sys.time()
  fit.lasso <- cv.glmnet(cbind(t, x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  fit.lasso.treat <- cv.glmnet(x=as.matrix(x), y=t, intercept=TRUE, family="binomial")
  
  activeX <- which(coef(fit.lasso.treat, s='lambda.1se')[-1] != 0)
  activeY <- which(coef(fit.lasso, s='lambda.1se')[-c(1,2)] != 0)
  DoubleActive = union(activeX, activeY)
  
  if (length(DoubleActive) == 0) {
    ## just do difference in means if no variables are included. No need for TMLE
    ModelDouble = lm(y ~ t)
    est = ModelDouble$coef[2]
    Interval = c(est - 1.96*sqrt(vcov(ModelDouble)[2,2]),
                            est + 1.96*sqrt(vcov(ModelDouble)[2,2]))
  } else {
    result1 <- tmle(Y=y,A=t,W=x[,DoubleActive],
                    g.SL.library = c("SL.glm"),
                    Q.SL.library = c("SL.glm"))
    
    est = result1$estimates$ATE$psi
    Interval = c(result1$estimates$ATE$CI[1],
                 result1$estimates$ATE$CI[2])
  }
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

DRmatch = function(y, t, x) {
  time1 = Sys.time()
  n = length(y)
  
  fit.x <- cv.glmnet(x, t, intercept = TRUE, family = 'binomial')
  fit.y <- cv.glmnet(x[t==0,], y[t==0], intercept = TRUE)
  fit.xy <- cv.glmnet(cbind(t,x), y, intercept = TRUE, penalty.factor = rep(c(0,1), c(1, p)))
  
  xh <- expit(predict(fit.x, newx = x))
  xhlogit <- predict(fit.x, newx = x)
  yh <- predict(fit.y, newx = x)
  
  yh = yh + rnorm(n, sd=0.00001)
  xhlogit = xhlogit + rnorm(n, sd=0.00001)
  
  fit.match.two = Match(Y=y, Tr=t, X=cbind(xhlogit,yh), estimand = "ATE", caliper=0.5)
  
  yh.xy = predict(fit.xy, newx=cbind(t,x))
  sigmaY2est = mean((yh.xy - y)^2)
  
  nTreatSig = table(fit.match.two$index.treated)^2
  varTreat = sum(nTreatSig)*sigmaY2est / length(fit.match.two$index.treated)^2
  nTreatCon = table(fit.match.two$index.control)^2
  varCon = sum(nTreatCon)*sigmaY2est / length(fit.match.two$index.control)^2
  varTot = varTreat + varCon
  seTot3 = sqrt(varTot)
  
  est = fit.match.two$est
  Interval = c(est - 1.96*seTot3,
               est + 1.96*seTot3)
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

DRB = function(y, t, x) {
  time1 = Sys.time()
  estLinear = DRbayes(y=y, t=t, x=x, nScans=5000, nBurn=3000, thin=2)
  
  est = estLinear$TreatEffect
  Interval = c(estLinear$TreatEffectCI[1],
               estLinear$TreatEffectCI[2])
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

HDC = function(y, t, x) {
  time1 = Sys.time()
  sslEB = SSL(y=y, z=t, x=x, nScans=5000, burn=3000, thin=2, lambda0="EB")
  
  est = sslEB$TreatEffect
  Interval = c(sslEB$TreatEffectCI[1],
               sslEB$TreatEffectCI[2])
  
  time2 = Sys.time()
  
  l = list(est = est,
           Interval = Interval,
           time = as.numeric(difftime(time2, time1, units="secs")))
  
  return(l)
}

expit = function(x) {
  return(exp(x) / (1 + exp(x)))
}
