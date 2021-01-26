# Compare geem2 with geem 
devtools::load_all()

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


if (FALSE) {
devtools::load_all()
i <- 1

if (res$waves[i] == "none") {
  theseWaves <- NULL 
} else {
  theseWaves <- resp[, res$waves[i]]
}

theseWeights <- resp[, res$weights[i]]
thisTreat <- res$treat[i]
thisCorStr <- res$corstr[i]

thisFormula <- update(useFormula, as.formula(paste(". ~ . + ", thisTreat)))

mm <- geeM2::geem2(thisFormula, data = resp, 
                   id = with(resp, interaction(center, id)),
                   family = "binomial", corstr = thisCorStr, Mv = thisMv,
                   corr.mat = thisCorr, waves = theseWaves,
                   weights = theseWeights)
debugonce(geem2)
}





if (FALSE) {

  debugonce(geem2)
mm_gm2 <- geem2(outcome ~ baseline + center + sex + age + I(age^2) + treat, #wmiss1, 
                 data = resp, 
                 id = with(resp, interaction(center, id)),
                     family = "binomial", corstr = "exchangeable",
                output = "geeglm")#,
                   #  waves = theseWaves,
#                     weights = weights1)

mm_gp <- geeglm(outcome ~ baseline + center + sex + age + I(age^2) + treat, #wmiss1,
                data = resp[, #complete.cases(resp$treatwmiss1), 
                            c("outcome", "baseline","sex", "age", "treat", #wmiss1", 
                              "id", "center")],#, "weights1")], 
              id = interaction(center, id),
              family = "binomial", corstr = "exchangeable")#, 
          #    waves = theseWaves[complete.cases(resp)],
             # weights = weights1)

mm_gp$model
mm_gp$contrasts
mm_gp$xlevels

nstop <- 0
for (i in 1:nstop) {
  print(i)
}


head(data.frame(mm_gp$weights, mm_gp$prior.weights))

nrow(resp[complete.cases(resp) * theseWeights != 0,]) - mm_gp$rank

mm_gp$df.residual

complete.cases(resp) * theseWeights != 0


mm_glm <- glm(outcome ~ baseline + center + sex + age + I(age^2) + treat, 
              data = resp, family = "binomial")

head(mm_glm$data)


head(mm_glm$model)

rr <- resp[,  c("outcome", "treat", "id", "center")]
neworder <- order(interaction(rr$id))



oldorder <- c(1:444)[neworder]
rr_new <- rr[neworder,]

head(rr_new, 20)
head(rr_new[order(neworder),], 20)


head(rr_new[1:44,])
neworder


td <- data.frame(originalid = 1:20, id = rep(letters[1:4], 5))
td

tdord <- order(td$id)
tdord
td2 <- td[tdord,]
tdord[td$originalid, ]
}
