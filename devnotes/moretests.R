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

m2 <- geeglm(outcome ~ treat + navarused, 
             data = resp2,
             id = interaction(center, id),
             family = "binomial", corstr = "exchangeable")

coef(m1); coef(m2)

summary(m1); summary(m2)
  
 # Conclusion: 
 # seems fine: missing information in "other" variables is ignored.

####################################################################################

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
data("respiratory")

#model
m <- geem2(outcome ~ treat, 
             data = respiratory,
             id = with(respiratory, interaction(center, id)),
             family = "binomial", corstr = "exchangeable",
           output = "geeglm")
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


####################################################################################
# Investigate handling of non-equidistant time points via waves argument
####################################################################################
set.seed(123)
exdat <- data.frame(x = rnorm(20))
exdat$y <- exdat$x + rnorm(20) 
exdat$id <- rep(1:5, each = 4)
exdat$time <- c(1, 2, 3, 4, 1, 2, 5, 6, 1, 2, 4, 3, 5, 6, 1, 10, 1, 2, 3, 4)
exdat <- exdat[sample(1:20, 20),] #scramble order

#corstr: exchangeable
m_exch_waves <- geem2(y ~ x, id = id, waves = time,
                      data = exdat, corstr = "exchangeable",
                      output = "geeglm")

m_exch_nowaves <- geem2(y ~ x, id = id,
                      data = exdat, corstr = "exchangeable",
                      output = "geeglm")

identicalLists(m_exch_waves, m_exch_nowaves, ignore = c("call"), chatty = TRUE) #OK

#corstr ar1:
m_ar1_waves_a <- geem2(y ~ x, id = id, waves = time,
                      data = exdat, corstr = "ar1",
                      output = "geeglm", testARG = 0)

m_ar1_waves_b <- geem2(y ~ x, id = id, waves = time,
                     data = exdat, corstr = "ar1",
                     output = "geeglm", testARG = 5)

m_ar1_waves_c <- geem2(y ~ x, id = id, waves = time,
                     data = exdat, corstr = "ar1",
                     output = "geeglm", testARG = 100)

m_ar1_waves_d <- geem2(y ~ x, id = id, waves = time,
                       data = exdat, corstr = "ar1",
                       output = "geeglm", testARG = NULL)

m_ar1_nowaves <- geem2(y ~ x, id = id,
                        data = exdat, corstr = "ar1",
                        output = "geeglm")

m_ar1_waves_a
m_ar1_waves_b
m_ar1_waves_c
m_ar1_waves_d
m_ar1_nowaves


m_ar1_waves
m_ar1_nowaves
 #different alpha estimates but similar beta estimates - promising.

debugonce(geem2)
m_ar1_waves

identicalLists(m_ar1_waves, m_ar1_nowaves, ignore = c("call"), chatty = TRUE) #OK
               