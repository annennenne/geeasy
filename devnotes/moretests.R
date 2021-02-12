library(geepack)
library(devtools)
load_all()
source("devnotes/helpers.R") 


####################################################################################
# Investigate output format when there is missing information in input 
# - should correspond to lm/glm standards
####################################################################################

set.seed(123)
exdat <- data.frame(x = rnorm(20))
exdat$y <- exdat$x + rnorm(20) 
exdat$x[sample(1:20, 5)] <- NA
exdat$id <- rep(1:5, each = 4)
exdat <- exdat[sample(1:20, 20),] #scramble order

lm0 <- lm(y ~ x, exdat)
glm0 <- glm(y ~ x, data = exdat)
gm20 <- geelm(y ~ x, id = id, data = exdat, corstr = "independence")
lm0; glm0; gm20 #same results

#lm and glm have dropped NAs in output - order is not preserved
length(lm0$residuals)
length(glm0$residuals)
length(lm0$fitted.values)
length(glm0$fitted.values)
length(gm20$fitted.values)

cbind(fitted(lm0), fitted(glm0), fitted(gm20))
 #identical outputs, same ordering


####################################################################################
# Test methods for geeM2s geepack output
####################################################################################

library(geepack)
library(geeasy)
data("respiratory")
respiratory$useid <- with(respiratory, interaction(center, id))

#model
m0 <- geelm(outcome ~ 1, 
           data = respiratory,
           id = useid,
           family = "binomial", corstr = "exchangeable")
m <- geelm(outcome ~ treat, 
             data = respiratory,
             id = useid,
             family = "binomial", corstr = "exchangeable")
m2 <- geelm(outcome ~ treat + sex + age + baseline, 
            data = respiratory,
            id = useid,
            family = "binomial", corstr = "exchangeable")
m2_ar1 <- geelm(outcome ~ treat + sex + age + baseline, 
            data = respiratory,
            id = useid,
            family = "binomial", corstr = "ar1")
m2_indep <- geelm(outcome ~ treat + sex + age + baseline, 
                data = respiratory,
                id = useid,
                family = "binomial", corstr = "independence")
mgp0 <- geeglm(outcome ~ 1, 
              data = respiratory,
              id = useid,
              family = "binomial", corstr = "exchangeable")
mgp <- geeglm(outcome ~ treat, 
              data = respiratory,
              id = useid,
              family = "binomial", corstr = "exchangeable")
mgp2 <- geeglm(outcome ~ treat + sex + age + baseline, 
              data = respiratory,
              id = useid,
              family = "binomial", corstr = "exchangeable")

#print
m

#summary & summary.print
summary(m)

#predict
predict(m)
predict(m, newdata = data.frame(treat = c("A", "P")))

#confint (NEW!)
confint(m)

#drop1
 drop1(m2)
 
#plot 
  #plot one model (geelm) 
  plot(m)
  
  #plot two models (geelm)
  plot(m, m2)
  
  #plot geelm and geeglm together
  plot(geelm = m2, geeglm = mgp2)
  
  #plot geelms with different corstr together for comparison
  plot(Exchangeable = m2, Independent = m2_indep,
       AR1 = m2_ar1)
  
  #test: more models than colors in palette 
  plot(a = m2, b = m2, c = m2, d = m2, 
       e = m2, f = m2, g = m2, h = m2, i = m2,
       j = m2, k =m2)
    #works

#QIC
QIC(m)

#getGEE
  getGEE(m, "beta")
  getGEE(m, "beta.se")
  getGEE(m, "alpha")
  getGEE(m, "nclusters")
  getME(m, "beta")
  getME(mgp, "beta") #works on geepack::geeglm output as well

#anova
  #zero covariates - doesn't work! 
  # - but doesn't work for original geeglm either (bug reported)
  anova(m0) 
  anova(mgp0)
  
  #one covariate - works
  anova(m)
  
  #multiple covariates - works
  anova(m2) 
  anova(mgp2)
  
  #comaparison of two nested models - works
  anova(m, m2) 
  anova(mgp, mgp2)




####################################################################################
# Handling of non-equidistant time points via waves argument
####################################################################################

set.seed(123)
exdat <- data.frame(x = rnorm(20))
exdat$y <- exdat$x + rnorm(20) 
exdat$id <- rep(1:5, each = 4)
exdat$time_nonequi <- c(1, 2, 3, 4, 1, 2, 5, 6, 1, 2, 4, 3, 5, 6, 1, 10, 1, 2, 3, 4)
exdat$time_equi <- c(1:4, 1:4,  1, 2, 4, 3,   2, 3, 1, 4,   1:4)
exdat <- exdat[sample(1:20, 20),] #scramble order

m_ar1_waves_nonequidist <- geelm(y ~ x, id = id, waves = time_nonequi,
                                 data = exdat, corstr = "ar1",
                                 output = "geeglm")

m_ar1_waves_equidist <- geelm(y ~ x, id = id, waves = time_equi,
                              data = exdat, corstr = "ar1",
                              output = "geeglm")

m_ar1_nowaves <-  geelm(y ~ x, id = id,
                        data = exdat, corstr = "ar1",
                        output = "geeglm")

#compare outputs
  m_ar1_waves_nonequidist
  m_ar1_waves_equidist
   #identical results (as expected)

  m_ar1_nowaves
    #different results than the other
    #this is also as expected, since the waves are not ordered 1:4
  
  
  
####################################################################################
## New functionality: use geepack as engine from geeasy::geelm
####################################################################################

set.seed(123)
exdat <- data.frame(x = rnorm(20))
exdat$y <- exdat$x + rnorm(20) 
exdat$id <- rep(1:5, each = 4)
exdat2 <- exdat[sample(1:20, 20),] #scramble order  
  
m_eng_geeasy_ordered <-  geelm(y ~ x, id = id,
                         data = exdat, corstr = "exchangeable")
m_eng_geepack_ordered <-  geelm(y ~ x, id = id,
                                data = exdat, corstr = "exchangeable",
                                engine = "geepack")
m_eng_geeasy_scrambled <-  geelm(y ~ x, id = id,
                               data = exdat2, corstr = "exchangeable")
m_eng_geepack_scrambled <-  geelm(y~ x, id = id,
                                data = exdat2, corstr = "exchangeable",
                                engine = "geepack")

# Compare geeasy engined results with/without scrambling
m_eng_geeasy_ordered
m_eng_geeasy_scrambled
  #Success: Identical results

# Compare geeasy and geepack engines (ordered data)
m_eng_geeasy_ordered
m_eng_geepack_ordered
  #Success: Similar results 

# Compare geeasy and geepack engines (scrambled data)
m_eng_geeasy_scrambled
m_eng_geepack_scrambled
  #Success: Similar results

# Compare geepack engined results with/without scrambling
m_eng_geepack_ordered
m_eng_geepack_scrambled
  #Success: Identical results



####################################################################################
## Test that geelm finds weights/waves/id in data, in global env, as expressions
####################################################################################

set.seed(123)
exdat <- data.frame(x = rnorm(20))
exdat$y <- exdat$x + rnorm(20) 
exdat$thisid <- rep(1:5, each = 4)
exdat$ws <- rep(seq(0.1, 1, 0.1), 2)
exdat$time <- rep(1:4, 5)

global_id <- exdat$thisid
global_ws <- exdat$ws
global_time <- exdat$time

# Test: Model with id/weights/waves given as names from data (unquoted)
geelm(y ~ x, data = exdat, 
      id = thisid, 
      weights = ws, 
      waves = time,
      corstr = "exchangeable")
  #works!

# Test: Model with id/weights/waves given as names from data (quoted)
geelm(y ~ x, data = exdat, 
      id = "thisid", 
      weights = "ws", 
      waves = "time",
      corstr = "exchangeable")
  #works!

# Test: Model with id/weights/waves given as variables from global env
geelm(y ~ x, data = exdat, 
      id = global_id, 
      weights = global_ws, 
      waves = global_time,
      corstr = "exchangeable")
  #works!

# Test: Model with id/weights/waves given as expressions using names from data
geelm(y ~ x, data = exdat, 
      id = thisid + 2, 
      weights = ws + 2, 
      waves = time + 2, 
      corstr = "exchangeable")
  #works!

# Test: Model with id/weights/waves given as expressions without names from data
geelm(y ~ x, data = exdat,
      id = rep(1:5, each = 4),
      weights = rep(seq(0.1, 1, 0.1), 2),
      waves = rep(1:4, 5),
      corstr = "exchangeable")
 #works!








