# geeasy
[![CRAN\_Release\_Badge](http://www.r-pkg.org/badges/version-ago/geeasy)](https://CRAN.R-project.org/package=geeasy)
![Download counter](http://cranlogs.r-pkg.org/badges/grand-total/geeasy)

R package with tools for fitting generalized linear models with
clustered observations using generalized estimating equations.

## Installation

geeasy is available on CRAN and can be installed as follows:

```{r}
install.packages("geeasy")
```

To install the development version of `geeasy` run the following
commands from within R (requires that the `devtools` package is already
installed)

```{r}
devtools::install_github("annennenne/geeasy")
```


## Fitting GEE models

`geeasy` fits generalized linear models on data with
correlated/clustered observations by use of generalized estimating
equations:

```{r}
library(geeasy)

# load data
data("respiratory")
respiratory$useid <- interaction(respiratory$center, respiratory$id)

# fit model
m <- geelm(outcome ~ treat + sex + age + baseline, data = respiratory,
           id = useid, family = "binomial", corstr = "exchangeable")
```

The syntax is similar to `glm()`, but a few additional arguments need to be specified:

* `id`: ID for identifying clusters. All observations with the same id
  are considered to belong to the same cluster.

* `corstr`: The correlation structure that is used within each
  cluster. Options include "independence" (the default, corresponding
  to no clustering), "exchangeable" (identical pair-wise correlations
  between all observations within a cluster) and more, see the
  documentation of `geelm()` for more details.


## Tools for working with GEE models

The package includes a selection of functions that can be used to
inspect and work with GEE models. These functions can be used
both with the output from `geelm()` *and* with the output of
`geeglm()` from the `geepack` R package.

The following functions are implemented in `geeasy`:

* `getGEE()`
* `plot()`
* `confint()`
* `drop1()`

A few more details about the two non-standard functions, `getGEE()`
and `plot()` are provided below. Furthermore, `geeasy` imports the
following functions from `geepack` that are also available:

* `summary()`
* `print()`
* `anova()`
* `QIC()`


**`getGEE()`**: 
```{r}
# Get parameter estimates:
getGEE(m, "beta")

# Get standard errors for parameter estimates: 
getGEE(m, "beta.se")

# Get estimated alpha (correlation structure parameter):
getGEE(m, "alpha")
```

This function was built to resemble the `getME()` function from
`lme4`. Note that it can also be accessed by calling `getME()`.

**`plot()`**:
```{r}
# Plot estimates and 95% confidence intervals for one geelm model
plot(m)

# Fit a new geelm model with AR1 correlation structure AND a glm 
# (corresponding to independent correlation structure)
m_ar1 <- geelm(outcome ~ treat + sex + age + baseline, 
               data = respiratory, id = useid,
               family = "binomial", corstr = "ar1")
m_glm <- glm(outcome ~ treat + sex + age + baseline, 
               data = respiratory, family = "binomial")
               
# Plot all three models together for easy comparison
plot(m, m_ar1, m_glm)
```

Note that this plotting function can also be accessed by calling
`plotEst()` and that this function allows for any number of models to
be plotted together, and it supports the model types `lm`, `glm`,
`geelm`, `geeglm`, `mice` and more.

## More options for `geelm()`

**Changing the output object:** `geelm()` can output a `geem` object,
resembling the output of `geem()` from the `geeM` package:

```{r}
m_outout_geem <- geelm(outcome ~ treat + sex + age + baseline, 
                       data = respiratory, id = useid,
                       family = "binomial", corstr = "exchangeable",
                       output = "geem")
```
This does not change the computations performed, only the output object. This means that the output will generally *not* be identical to that of `geeM::geem()`.

**Changing the estimation engine:** `geelm()` allows for choosing to use `geepack` as its computational engine as follows:

```{r}
m_engine_geepack <- geelm(outcome ~ treat + sex + age + baseline, 
                   data = respiratory, id = useid,
                   family = "binomial", corstr = "exchangeable",
                   engine = "geepack")
```
Note that this does **not** mean that the id variable is handled as in `geepack`: Clusters are still constructed by assigning observations with identical values of `id` to the same cluster.


## Credit

The `geeasy` package is based on a modified version of the `geeM` package and the main estimation code was hence written by Lee McDaniel and Nick Henderson. 

The package was modified, updated and extended by Anne Helby Petersen. 

Claus Ekstrøm has contributed additional code. 

Søren Højsgaard is maintainer of the `geeasy` package. 


## Bugs & requests

If you find bugs or have a request for a new feature, please [open an issue](https://github.com/annennenne/geeasy/issues).
