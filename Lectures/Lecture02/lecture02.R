## ----setup, include=FALSE, cache=FALSE----------------------------------------
#Load packages and initialize options
opts_chunk$set(fig.path='./knitr-figs/lecture01-',echo=FALSE, cache=TRUE,results='hide', warning=FALSE, fig.width=8,fig.height=4.5,size='scriptsize',error=FALSE,tidy=FALSE)
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
###  Rnw File for lecture 02
###
### History:
###  -- 2020-06-10 file created as part of SmitSjuk
###  -- 2020-02-27 modified as part of STA427
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(RColorBrewer)


## ----echo=TRUE,results='markup',tidy=FALSE------------------------------------
######################################################################
# Likelihood function for the Reed-Frost model
#
# Parameters:
#  w.logit - logit(w) to have unrestricted parameter space
#  x       - vector containing the number of susceptibles at each time
#  y       - vector containing the number of infectious   at each time
#
######################################################################

l <- function(w.logit,x,y) {
  if (length(x) != length(y)) { stop("x and y need to be the same length") }
  K <- length(x)
  w <- plogis(w.logit)
  p <- 1 - (1-w)^y
  return(sum(dbinom( y[-1], size=x[-K], prob=p[-K],log=TRUE)))
}

# Epidemic D in Table 4.1 of Daley and Gani (1999), assuming all susceptibles got infected
y <- c(1, 4, 14, 10, 1, 0)
x <- numeric(length(y))
x[1] <- sum(y[-1])
x[2:length(x)] <- x[1]-cumsum(y[2:length(y)])

mle <- optim(par=0,fn=l,method="BFGS",x=x,y=y,control=list(fnscale=-1),hessian=TRUE)
# Maximum likelihood estimator
(w.hat <- plogis(mle$par))

## ----echo=FALSE, results="hide"-----------------------------------------------
# 95% confidence interval
(w.95ci <- plogis( mle$par + c(-1,1)*qnorm(0.975)* sqrt(-1/as.vector(mle$hess))))


## ----fig.align="center"-------------------------------------------------------
x0 <- seq(sum(y[-1]), sum(y[-1])+20, 1)
w.hat <- numeric(length(x0))
log.lik <- numeric(length(x0))
for(i in 1:length(x0)){
    x[1] <- x0[i]
    x[2:length(x)] <- x[1]-cumsum(y[-length(y)])
    
    mle <- optim(par=0,fn=l,method="BFGS",x=x,y=y,control=list(fnscale=-1),hessian=TRUE)
    w.hat[i] <- plogis(mle$par)
    log.lik[i] <- mle$value
}

par(mfrow=c(1,1))
plot(x0, log.lik, type="b", ylab="profile log likelihood", axes=TRUE)


## ----readcsfvoutbreak---------------------------------------------------------
csfv <- read.table(file.path("..","Lecture01","Data","csfv.txt"),col.names=c("t","I"))
pal <- brewer.pal(3,"Set1")


## ----echo=FALSE, results="hide"-----------------------------------------------
# Load pkg.
suppressPackageStartupMessages(library(deSolve))

##############################################################################
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list containing dS(t)/dt and dI(t)/dt
##############################################################################

sir <- function(t,y, parms) {
  beta <- parms[1]
  gamma <- parms[2]
  S <- y[1]
  I <- y[2]
  return(list(c(S=-beta*S*I,I=beta*S*I-gamma*I)))
}

## ----echo=FALSE---------------------------------------------------------------
N <- 21500
#sumy <- sum(csfv$I)
sumy <- 429
f <- sumy/N
R0 <- -log(1-f)/f


## ----echo=TRUE----------------------------------------------------------------
######################################################################
#Least-squares fit
######################################################################

ll.gauss <- function(theta, take.sqrt=FALSE) {
  #Solve ODE using the parameter vector theta
  res <- lsoda(y=c(N-1,1), times=csfv$t, func=sir, parms=exp(theta))
  #Squared difference?
  if (take.sqrt==FALSE) {
    return(sum(dnorm(csfv$I,mean=res[,3],sd=1,log=TRUE)))
  } else { 
    return(sum(dnorm(sqrt(csfv$I),mean=sqrt(abs(res[,3])),sd=1,log=TRUE)))
  }
}


## ----echo=TRUE----------------------------------------------------------------
#Determine MLE
N <- 21500
mle <- optim(log(c(0.00002,3)), fn=ll.gauss,control=list(fnscale=-1))

#Show estimates and resulting R0 estimate
beta.hat <- exp(mle$par)[1]
gamma.hat <- exp(mle$par)[2]
R0.hat <- beta.hat*N/gamma.hat


## ----echo=TRUE, results="show"------------------------------------------------
mu <- lsoda(y=c(N-1,1), times=csfv$t, func=sir,parms=exp(mle$par))
head(mu, n=3)


## ----echo=FALSE---------------------------------------------------------------
mle2 <- optim(log(c(0.00002,3)), fn=ll.gauss, take.sqrt=TRUE, control=list(fnscale=-1))
beta.hat2 <- exp(mle2$par)[1]
gamma.hat2 <- exp(mle2$par)[2]
R0.hat2 <- beta.hat2*N/gamma.hat2
mu2 <- lsoda(y=c(N,1), times=csfv$t, func=sir,parms=exp(mle2$par))


## ----fig.width=7,fig.height=4.0, fig.align="center"---------------------------
pal <- brewer.pal(4,"Set1")
matplot(mu[,1],cbind(csfv$I,mu[,3],mu2[,3]),type="l",lwd=3,lty=1,ylab="No. infectious herds",xlab="time (weeks after first infection)",col=pal)
legend(x="topright",c("CSFV outbreak","LS fit", "LS-sqrt fit"), lty=1,col=pal,lwd=3)


## ----echo=FALSE, warning=FALSE------------------------------------------------
######################################################################
# Poisson likelihood
######################################################################

ll.pois <- function(theta) {
  #Solve ODE using the parameter vector theta
  res <- lsoda(y=c(N-1,1), times=csfv$t, func=sir, parms=exp(theta))
  #Poisson likelihood
  return(sum(dpois(round(csfv$I),lambda=res[,3],log=TRUE)))
}

#Determine MLE
mle <- optim(mle$par, fn=ll.pois,control=list(fnscale=-1))

#Show results
(beta.hat <- exp(mle$par)[1])
(gamma.hat <- exp(mle$par)[2])
(R0.hat <- beta.hat*N/gamma.hat)
muP <- lsoda(y=c(N-1,1), times=csfv$t, func=sir,parms=exp(mle$par))

## ----fig.width=7,fig.height=4.0-----------------------------------------------
matplot(mu[,1],cbind(csfv$I,mu[,3],mu2[,3], muP[,3]),type="l",lwd=3,lty=1,ylab="No. infectious herds",xlab="time (weeks after first infection)",col=pal)
legend(x="topright",c("CSFV outbreak","LS fit", "LS-sqrt fit", "Poisson fit"),lty=1,col=pal,lwd=3)


## ----echo=FALSE,eval=FALSE----------------------------------------------------
## #Estimate R0 from final size data
## f <- 429/N
## (R0.hat.finalsize <- -log(1-f)/f)

