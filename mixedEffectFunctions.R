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
  
  y <- intercept + randomEffectsItem[, "intercept"][id] + randomEffectsCondition[, "intercept"][id] +
    (fixedEffectCondition + randomEffectsCondition[, "slope"][id]) * fixedEffectConditionVector[cond] +
    (fixedEffectItem + randomEffectsItem[, "slope"][id]) * fixedEffectItemVector[item] + err
  
  mydat <- data.frame(y = y,
                      cond = as.factor(cond),
                      item = as.factor(item),
                      id = as.factor(id))
  return(mydat)
}
