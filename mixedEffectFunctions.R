library(MASS)
library(lme4)
library(BayesFactor)
library(afex)
library(lmerTest)
# library(BFEffects)

myDataSim <- function(fixedEffectItem = 0, randomSlopeItemSigma = 0, randomInterceptItemSigma = 0, randomRhoSlopeInterceptItem = 0,
                      fixedEffectCondition = 0, randomSlopeConditionSigma = 0, randomInterceptConditionSigma = 0, randomRhoSlopeInterceptCondition = 0,
                      errorSigma = 1, errorSigBetween = 0, nSub = 30, nItem = 2, nCond = 2) {
  
  I <- nSub # n per item per condition
  id <- 1:I
  intercept <- 0
  
  M <- nItem # items
  J <- nCond # conditions
  grandN <- I * J * M
  cond <- rep(1:J, each = M*I)
  item <- rep(1:M, J*I)
  id <- rep(rep(1:I, each = M), J)
  fixedEffectConditionVector <- seq(-1/nCond, 1/nCond, length.out = nCond) 
  fixedEffectItemVector <- seq(-1/nItem, 1/nItem, length.out = nItem) 
  
  # setup random effects  
  randomEffectMatItem <- diag(c(randomInterceptItemSigma, randomSlopeItemSigma))
  randomEffectMatItem[1, 2] <- randomEffectMatItem[2, 1] <- randomRhoSlopeInterceptItem * sqrt(randomInterceptItemSigma) * sqrt(randomSlopeItemSigma)
  
  randomEffectMatCondition <- diag(c(randomInterceptConditionSigma, randomSlopeConditionSigma))
  randomEffectMatCondition[1, 2] <- randomEffectMatCondition[2, 1] <- randomRhoSlopeInterceptCondition * sqrt(randomInterceptConditionSigma) * sqrt(randomSlopeConditionSigma)
  
  err <- rnorm(grandN, 0, errorSigma)

  randomEffectsItem <- mvrnorm(n = I, mu = c(0,0) , Sigma = randomEffectMatItem)
  randomEffectsCondition <- mvrnorm(n = I, mu = c(0,0) , Sigma = randomEffectMatCondition)
  colnames(randomEffectsItem) <- colnames(randomEffectsCondition) <- c("intercept", "slope")
  
  y <- intercept + randomEffectsItem[,"intercept"][id] + randomEffectsCondition[,"intercept"][id] +
    (fixedEffectCondition + randomEffectsCondition[,"slope"][id]) * fixedEffectConditionVector[cond] +
    (fixedEffectItem + randomEffectsItem[,"slope"][id]) * fixedEffectItemVector[item] + err
  
  mydat <- data.frame(y = y,
                      cond = as.factor(cond),
                      item = as.factor(item),
                      id = as.factor(id))
  return(mydat)
}


compute4BF=function(dat){
  y=dat$y
  n=length(y)
  X0=rep(1,n)
  I <- length(unique(dat$id))
  Xsub=matrix(0,nrow=n,ncol=I)
  for (i in 1:n) Xsub[i, dat$id[i]]=1
  X1 = sign(as.numeric(dat$cond) - 1.5) * (1/sqrt(2))
  XsubCond=matrix(0,nrow=n,ncol=I)
  for (i in 1:n) XsubCond[i, dat$id[i]] <- sign(as.numeric(dat$cond[i]) - 3/2) * (1/sqrt(2))
  
  bf=1:4
  #Mod 1 Strong Null, no condition effects
  X=cbind(X0,Xsub)
  g=rep(0,I)
  bf[1]=nWayAOV(y,X,g,rscale=1)$bf
  #Mod 2 No subject-by-treatment interactions
  X=cbind(X0,Xsub,X1)
  g=c(rep(0,I),1)
  bf[2]=nWayAOV(y,X,g,rscale=c(1,1))$bf
  #Mod 3  Zero-Centered subject-by-treatment interactions (silly model)
  X=cbind(X0,Xsub,XsubCond)
  g=c(rep(0,I),rep(1,I))
  bf[3]=nWayAOV(y,X,g,rscale=c(1, 1))$bf
  #Mod 4 Main effect with subject-by-treatment interactions
  X=cbind(X0,Xsub,X1,XsubCond)
  g=c(rep(0,I),1,rep(2,I))
  bf[4]=nWayAOV(y,X,g,rscale=c(1, 1, 1))$bf
  out=bf
  names(out)=c("M0","M2", "M3", "M4")
  return(out)
}  