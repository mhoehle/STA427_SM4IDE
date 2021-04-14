## ----setup, include=FALSE-----------------------------------------------------
#Load packages and initialize options
opts_chunk$set(fig.path='./knitr-figs/lecture06-',echo=FALSE, cache=FALSE,results='hide', warning=FALSE, fig.width=8,fig.height=4.5,size='tiny',error=FALSE,tidy=FALSE)
options(replace.assign=TRUE,width=80)


## ----echo=FALSE,results='hide',message=FALSE----------------------------------
######################################################################
### Author: Michael Höhle <http://www.math.su.se/~hoehle>
### Course: Statistical Methods in Infectious Disease Epidemiology
###         (STA427 - Spring 2021) at the
###         Department of Biostatistics,
###         Epidemiology, Biostatistics and Prevention Institute,
###         University of Zurich, Switzerland
###
### License:
### Slides content is licensed under a <a rel="license"
### href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
### Attribution-ShareAlike 4.0 International License</a>.
### The code content is available under a <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a>
### license.
###
### Description:
###  Rnw File for lecture 07
###
### History:
###  -- 2021-04-10 file created
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(lubridate)
library(stringr)
library(RColorBrewer)
library(surveillance)


## ----make_carnival, echo=FALSE, results="markup"------------------------------
carnival <- bind_rows(
    data.frame(infected=0, n= 4, cs_source="asymptomatic"),
    data.frame(infected=4, n=28, cs_source="symptomatic")
  ) %>%
  mutate(across(cs_source, as.factor)) %>%
  mutate(cs_source = relevel(cs_source, ref="asymptomatic"))
carnival

## ----glm_carnival, echo=TRUE, results="markup"--------------------------------
m_glm <- glm( cbind(infected, n-infected) ~ 1 + cs_source, data=carnival, family=binomial )
confint(m_glm)


## ----elrm_carnival, results="hide", message=FALSE, echo=TRUE, cache=TRUE------
#devtools::install_github(repo="https://github.com/cran/elrm.git")
m <- elrm::elrm( infected/n ~ cs_source, interest=~cs_source, dataset=carnival, r=2, iter=1.5e4, burnIn=5e3)

## ----elrm_carnival_results, echo=TRUE, results="markup"-----------------------
c(hat=as.numeric(exp(m$coeffs)), exp(m$coeffs.ci)) %>% unlist


## ----riskscore, echo=TRUE, results="markup", cache=TRUE-----------------------
PropCIs::riskscoreci(4, 28, 0, 4, conf.level=0.95)


## ----rrci_bayes, echo=TRUE, results="markup", cache=TRUE----------------------
PropCIs::rrci.bayes(4, 28, 0, 4, a=1/2, b=1/2, c=1/2, d=1/2, conf.level=0.95)

## ----eval=FALSE---------------------------------------------------------------
## #Manual version
## X <- rbeta(1e6, 1/2 + 4, 1/2 + 28-4)
## Y <- rbeta(1e6, 1/2 + 0, 1/2 + 4)
## Z <- X/Y
## quantile(Z, prob=c(0.025, 0.975))


## ----make_carnival_oc---------------------------------------------------------
carnival_oc <- bind_rows(
    data.frame(infected=0, n=22, cs_source="asymptomatic"),
    data.frame(infected=3, n=25, cs_source="symptomatic_phaseunknown_or_both"),
    data.frame(infected=15, n=72, cs_source="symptomatic_presymptomatic_only"),
    data.frame(infected=2, n=29, cs_source="symptomatic_symptomatic_only"),
  ) %>%
  mutate(across(cs_source, as.factor)) %>%
  mutate(cs_source = relevel(cs_source, ref="asymptomatic"))


## ----elrm_oc, cache=TRUE, results="hide", message=FALSE, warning=FALSE--------
m_oc <- elrm::elrm( infected/n ~ cs_source, interest=~cs_source, dataset=carnival_oc, r=2, iter=1.5e6, burnIn=1e4)


## ----elrm_oc_results, results="asis"------------------------------------------
tab <- cbind(hat=as.numeric(exp(m_oc$coeffs)[-1]), exp(m_oc$coeffs.ci))
tab <- tab %>% mutate(`Clinical symptoms source` = str_replace_all(names(tab), "cs_source", "")) %>%
  as_tibble() %>%
  select(`Clinical symptoms source`, everything()) %>%
  mutate(upper = str_replace(upper, "Inf", "$\\\\infty$"))
xtable::xtable(tab) %>% print(include.rownames=FALSE, sanitize.text.function=function(x) x)


## ----test_coarstedatatools, echo=TRUE, results="markup"-----------------------
# Simple dataset with 3 individuals
# type: 0 = doubly interval censored, 1=single interval censored, 2=exact
dat <- data.frame(EL = c(1,2,3), ER=c(2,3,3), SL=c(10,4,9), SR=c(12,7,9), type=c(0,0,2))

# Fit log-normal distribution to the data
coarseDataTools::dic.fit(dat=dat, dist="L")

