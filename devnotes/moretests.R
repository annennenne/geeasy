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

