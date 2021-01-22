# Compare geem2 with geem 

# Scenario 1: Respiratory data

library(geepack)
data(respiratory)


resp <- respiratory

# Add extra variables to easily test different scenarios
set.seed(1234)

# No weights (all weight 1)
resp$noweights <- rep(1, nrow(resp))

# 0/1 weights
resp$weights1 <- sample(c(0,1), nrow(resp), prob = c(0.05, 0.95), replace = TRUE)

# cont. weights (in (0,1))
resp$weights2 <- runif(nrow(resp))

# missing information scattered with no clusters being fully missing
set.seed(1222)
resp$treatwmiss1 <- resp$treat
resp$treatwmiss1[sample(c(TRUE, FALSE), nrow(resp), prob = c(0.1, 0.9), replace = TRUE)] <- NA

# missing information scattered + 1st cluster fully missing
resp$treatwmiss2 <- resp$treatwmiss1
resp$treatwmiss2[1:4] <- NA

# only first cluster fully missing
resp$treatwmiss3 <- resp$treat
resp$treatwmiss3[1:4] <- NA

# one-step waves variable (equidistant)
resp$timeequi <- rep(c(2, 4, 3, 1), nrow(resp)/4)

# waves variable with gaps
resp$timegaps <- rep(c(7, 4, 5, 1), nrow(resp)/4)


allCorStr <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", 
               "fixed", "userdefined")

thisMv <- 2 #for m-dep. 
thisCorr <- matrix(c(1, 2, 0, 0, #for fixed/userdefined
                     2, 1, 3, 0,
                     0, 3, 1, 4,
                     0, 0, 4, 1), byrow = TRUE, 
                   nrow = 4, ncol = 4)

res <- data.frame(corstr = rep(allCorStr, each = 3*4*3), 
                  weights = rep(c("noweights", "weights1", "weights2"), 7*4*3),
                  treat = rep(rep(c("treat", "treatwmiss1", "treatwmiss2", "treatwmiss3"), each = 3*3), 7),
                  waves =  rep(rep(c("none", "timeequi", "timegaps"), each = 3), 7*4),
                  identical = rep(NA, 7*3*3*4),
                  numdiff = rep(NA, 7*3*3*4),
                  error_original = rep(NA, 7*3*3*4),
                  error_new = rep(NA, 7*3*3*4))

source("devnotes/helpers.R") #load list comparison function


useFormula <- formula(outcome ~ baseline + center + sex + age + I(age^2))

for (i in 1:nrow(res)) {
  if (res$waves[i] == "none") {
    theseWaves <- NULL 
  } else {
    theseWaves <- resp[, res$waves[i]]
  }
  
  theseWeights <- resp[, res$weights[i]]
  thisTreat <- res$treat[i]
  thisCorStr <- res$corstr[i]
  
  thisFormula <- update(useFormula, as.formula(paste(". ~ . + ", thisTreat)))
  
  set.seed(123)
  mm_original <- tryCatch(expr = {geeM::geem(thisFormula, data = resp, 
                            id = with(resp, interaction(center, id)),
                            family = "binomial", corstr = thisCorStr, Mv = thisMv,
                            corr.mat = thisCorr, waves = theseWaves,
                            weights = theseWeights)}, error = function(e) 0)
  if (identical(mm_original, 0)) {
    res$error_original[i] <- TRUE
  } else {
    res$error_original[i] <- FALSE
  }

  set.seed(123)
  mm_dev <- tryCatch({geeM2::geem2(thisFormula, data = resp, 
                         id = with(resp, interaction(center, id)),
                         family = "binomial", corstr = thisCorStr, Mv = thisMv,
                         corr.mat = thisCorr, waves = theseWaves,
                         weights = theseWeights)}, error = function(e) 0)
  
  if (identical(mm_dev, 0)) {
    res$error_new[i] <- TRUE
  } else {
    res$error_new[i] <- FALSE
  }
  
  if (!res$error_original[i] & !res$error_new[i]) {
#    idres <- identicalLists(mm_original, mm_dev, 
#                            ignore = c("call", "FunList", "dropped"))
    idres <- identicalLists(mm_original, mm_dev, 
                            ignore = c("call", "FunList", "dropped",
                                       "clusz", "X", "offset", "eta", 
                                       "dropped", "weights", "terms", 
                                       "y", "formula"))
    if (is.logical(idres)) {
      res$identical[i] <- idres
    } else {
      res$identical[i] <- FALSE
      res$numdiff[i] <- idres
    }
  }
}


View(res)
