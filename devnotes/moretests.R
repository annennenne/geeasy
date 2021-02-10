library(geepack)
library(devtools)
load_all()
source("devnotes/helpers.R") 

####################################################################################
# Investigate behavior of geepack::geeglm when data contains NA (but not 
# in variables used for)
####################################################################################

data("respiratory")
resp2 <- respiratory
resp2$varnotused <- sample(c(0, 1, NA), nrow(respiratory),
                           replace = TRUE)
resp2$navarused <- c(rep(NA, 4), sample(c(0,1), nrow(respiratory) - 4,
                                        replace = TRUE))

m1 <- geeglm(outcome ~ treat, 
                data = respiratory,
                id = interaction(center, id),
                family = "binomial", corstr = "exchangeable")

m2 <- geeglm(outcome ~ treat, 
             data = resp2,
             id = interaction(center, id),
             family = "binomial", corstr = "exchangeable")

coef(m1); coef(m2)

summary(m1); summary(m2)
  
 # Conclusion: 
 # seems fine: missing information in "other" variables is ignored.


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
exdat$helpnum <- 1:20

lm0 <- lm(y ~ x + helpnum, exdat)
glm0 <- glm(y ~ x + helpnum, data = exdat)
gm20 <- geem2(y ~ x + helpnum, id = id, data = exdat, corstr = "independence",
             output = "geeglm")
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
library(geek)
data("respiratory")
respiratory$useid <- with(respiratory, interaction(center, id))

#model
m0 <- geelm(outcome ~ 1, 
           data = respiratory,
           id = useid,
           family = "binomial", corstr = "exchangeable",
           output = "geeglm")
m <- geelm(outcome ~ treat, 
             data = respiratory,
             id = useid,
             family = "binomial", corstr = "exchangeable",
           output = "geeglm")
m2 <- geelm(outcome ~ treat + sex + age + baseline, 
            data = respiratory,
            id = useid,
            family = "binomial", corstr = "exchangeable",
            output = "geeglm")
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
 library(MESS) 
 drop1(m2)
 
#add1
 #todo? 

#plot 
  #use plotEstimates?

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
  #zero covariates - doesn't work 
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