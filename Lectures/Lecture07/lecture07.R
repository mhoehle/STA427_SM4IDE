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


## ----echo=FALSE, results="markup"---------------------------------------------
carnival <- bind_rows(
    data.frame(infected=0, n= 4, cs_source="asymptomatic"),
    data.frame(infected=4, n=28, cs_source="symptomatic")
  ) %>%
  mutate(across(cs_source, as.factor)) %>%
  mutate(cs_source = relevel(cs_source, ref="asymptomatic"))
carnival

## ----echo=TRUE, results="markup"----------------------------------------------
m_glm <- glm( cbind(infected, n-infected) ~ 1 + cs_source, data=carnival, family=binomial )
confint(m_glm)


## ----results="hide", message=FALSE, echo=TRUE, cache=TRUE---------------------
#devtools::install_github(repo="https://github.com/cran/elrm.git")
m <- elrm::elrm( infected/n ~ cs_source, interest=~cs_source, dataset=carnival, r=2, iter=1.5e4, burnIn=5e3)

## ----echo=TRUE, results="markup"----------------------------------------------
c(hat=as.numeric(exp(m$coeffs)), exp(m$coeffs.ci)) %>% unlist

