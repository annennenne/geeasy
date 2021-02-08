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
             output = "geeglm", testARG =NULL)
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
data("respiratory")

#model
m <- geem2(outcome ~ treat, 
             data = respiratory,
             id = with(respiratory, interaction(center, id)),
             family = "binomial", corstr = "exchangeable",
           output = "geeglm", testARG = NULL)
mgp <- geeglm(outcome ~ treat, 
              data = respiratory,
              id = with(respiratory, interaction(center, id)),
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

#drop1/add1
 #find in MESS?

#plot 
  #use plotEstimates?

#QIC
QIC(m)

#getME 
  #to do

#anova
  #to do


####################################################################################
# Hndling of non-equidistant time points via waves argument
####################################################################################
set.seed(123)
exdat <- data.frame(x = rnorm(20))
exdat$y <- exdat$x + rnorm(20) 
exdat$id <- rep(1:5, each = 4)
exdat$time_nonequi <- c(1, 2, 3, 4, 1, 2, 5, 6, 1, 2, 4, 3, 5, 6, 1, 10, 1, 2, 3, 4)
exdat$time_equi <- c(1:4, 1:4,  1, 2, 4, 3,   2, 3, 1, 4,   1:4)
exdat <- exdat[sample(1:20, 20),] #scramble order

m_ar1_waves_nonequidist <- geem2(y ~ x, id = id, waves = time_nonequi,
                                 data = exdat, corstr = "ar1",
                                 output = "geeglm", testARG = NULL)

m_ar1_waves_equidist <- geem2(y ~ x, id = id, waves = time_equi,
                              data = exdat, corstr = "ar1",
                              output = "geeglm", testARG = NULL)

m_ar1_nowaves <-  geem2(y ~ x, id = id,
                        data = exdat, corstr = "ar1",
                        output = "geeglm", testARG = NULL)

#compare outputs
  m_ar1_waves_nonequidist
  m_ar1_waves_equidist
   #identical results (as expected)

  m_ar1_nowaves
    #different results than the other
    #this is also as expected, since the waves are not ordered 1:4