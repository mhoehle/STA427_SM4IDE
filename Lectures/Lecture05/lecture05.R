## ----setup, include=FALSE, cache=FALSE----------------------------------------
#Load packages and initialize options
opts_chunk$set(fig.path='./knitr-figs/lecture05-',echo=FALSE, cache=TRUE,results='hide', warning=FALSE, fig.width=8,fig.height=4.5,size='tiny',error=FALSE,tidy=FALSE)
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
###  Rnw File for lecture 05
###
### History:
###  -- 2021-03-10 modified as part of STA427
###  -- 2012-
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(lubridate)
library(stringr)
library(RColorBrewer)
library(surveillance)


## ----SALK---------------------------------------------------------------------
salk <- t(matrix(c(57,200745,142,201229),2,2))
dimnames(salk) <- list(c("vacc","non-vacc"),c("case","nocase"))


## ----echo=TRUE,results='verbatim'---------------------------------------------
salk <- t(matrix(c(57,200745,142,201229),2,2))
dimnames(salk) <- list(c("vacc","non-vacc"),c("case","nocase"))
require("Epi")
tab <- twoby2(salk,F.lim=1e9,print=FALSE)
(RR <- tab$measures["             Relative Risk:",])
(VE <- 1-RR[c(1,3,2)])


## ----echo=TRUE----------------------------------------------------------------
N <- function( VT, VL, pi, alpha=0.05, beta=0.2) {
  (qnorm(alpha/2) + qnorm(1-beta))^2 * (1-pi*VL)^2 * (1-VT) / ( pi*(1-pi)*(VT-VL)^2)
}


## ----echo=TRUE,results='verbatim'---------------------------------------------
VT<-seq(0.5,0.9,by=0.1)
t(sapply(c(0.05,0.1,0.15), function(delta) {
  structure(round(N( VT=VT,VL=VT-delta, pi=0.5,alpha=0.05,beta=0.5)),names=VT)
}))


## ----echo=TRUE,results='verbatim'---------------------------------------------
#Load data on measles notifications and vaccine coverage in Leeds
(measles <- read.csv2(file=file.path("Data", "farrington1993-table1.csv")))
#Create offset from PCV
measles$o <- with(measles, log(Coverage / (100-Coverage)))
#Create factor representing the birth cohort
measles$bc <- as.factor(measles$Birth.cohort)
m <- glm( cbind(Vaccinated, Cases-Vaccinated) ~ offset(o), family=binomial,data=measles)
m.bc <- glm( cbind(Vaccinated, Cases-Vaccinated) ~ offset(o) + bc, family=binomial,data=measles)
#Compare models
(p <- anova(m, m.bc,test="Chisq")$"Pr(>Chi)"[2])
#Vaccine coverage
(VE <- (1-exp(as.numeric(coef(m)))))
#Confidence intervals
sort(1-exp(confint(m)))

