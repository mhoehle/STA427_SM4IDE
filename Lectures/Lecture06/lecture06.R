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
###  Rnw File for lecture 06
###
### History:
###  -- 2021-03-24 modified as part of STA427
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(lubridate)
library(stringr)
library(RColorBrewer)
library(surveillance)


## -----------------------------------------------------------------------------
library("surveillance")
data(salmNewport)

#Load Salmonella Newport data
#load("../Data/sNewport.RData")
#Load extra plotting functionality
######################################################################
# Helper function to make a plot showing year and week of a surveillance
# time series in interpretable fashion.
######################################################################

plotOne <- function(sts,start,now,doSurv=TRUE,extraYlab="",ylim=c(0,max(observed(sts))),dy.alarm=-0.5) {
  #Remove superfluous margin
  par(mar=c(5,4,0.2,1))
  #Put those value not up to now to zero
  observed(sts)[epoch(sts) > now] <- 0

  plot(sts,dx.upperbound=0,legend.opts=NULL,ylab=iconv(paste("No. reported cases",extraYlab,sep=""),"utf-8","latin1"),main="",axes=FALSE,xlab="Year/Reporting Week",xlim=c(0,nrow(sts)),ylim=ylim)

  week <- isoWeekYear(epoch(sts))$ISOWeek
  year <- isoWeekYear(epoch(sts))$ISOYear
  #Where the by 10 divisible weeks
  is.subset1 <- week %% 10 == 0 | week == 1
  #Where the by 5 divisible weeks
  is.subset2 <- week %% 5 == 0
  #Where the last week of the year
  is.year <- week %% 52 == 0
  #Show axis with this
  axis(1,line=1,cex.axis=0.8,labels=NA,lwd.ticks=0)
  axis(1,(1:length(week))[is.subset1],label=week[is.subset1],line=1,cex.axis=0.8)
  axis(1,(1:length(week))[is.subset2],label=rep("",sum(is.subset2)),line=1,tck=-0.01)
  axis(1,(1:length(week))[is.year]+0.5,label=rep("",sum(is.year)),line=1,tcl=1.3,lwd.ticks=1)
  #Where to plot the years
  where <- tapply(1:length(year),factor(year),mean)
  axis(1,where,label=unique(format(epoch(sts),"%Y")),line=-1.4,lwd.ticks=NA,tick=FALSE)
  axis(2)

  #Show today
  idx <- which(epoch(sts) >= start & epoch(sts) <= now)
  #Upper bound
  points(max(idx),0,pch=20,cex=1.5,col="blue")

  #Do surveillance for this time point
  if (doSurv) {
    s.sts <- farrington(sts, control=list(range=idx,alpha=0.001,b=2,w=3,trend=FALSE,correct=FALSE,limit54=c(0,1)))
    for (i in 1:length(idx)) {
      lines( c(idx[i]-1,idx[i])+0.5, rep(upperbound(s.sts)[i,],2),lwd=3,col="red")
      if (i<length(idx)) {
        lines( rep(idx[i]+0.5,2), upperbound(s.sts)[c(i,i+1),],lwd=1,col="red")
      }
    }
    #Show alarms
    par(xpd=NA)
    points(idx[alarms(s.sts)>0],rep(dy.alarm,sum(alarms(s.sts)>0)),pch=24,cex=1,col="red",lwd=1)
  }

  invisible()
}

sNewport <- aggregate(salmNewport, by = "unit")
#Make a reduced time series for illustration
sNewportAll <- sNewport[epoch(sNewport) >= as.Date("2009-01-01") & epoch(sNewport) <= as.Date("2011-12-31"),]

#From where to start
start <- as.Date("2011-09-26")+7
surveillance::isoWeekYear(start)


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 0*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 1*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 2*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 3*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 4*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 5*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 5*7,doSurv=TRUE,ylim=c(0,20))
legend(x="topleft",expression(q[0.95]* " of predictive distribution"),col="red",lwd=2L,lty=1)


## ----echo=FALSE---------------------------------------------------------------
data(momo)


## ----echo=FALSE---------------------------------------------------------------
plot(momo[year(momo)>=2000,],type=observed ~ time | unit,par.list=list(mar=c(4,4,1,1)),ylab="No. deaths")


## ----echo=FALSE---------------------------------------------------------------
momo2 <- momo
momo2@observed <- observed(momo)/population(momo) * 100000
plot(momo2[year(momo2)>=2000,],ylab="Deaths per 100.000",type=observed ~ time | unit,par.list=list(mar=c(4,4,1,1)))


## ----PREDINT------------------------------------------------------------------
suppressPackageStartupMessages(library(gamlss))
set.seed(123)
n <- 5
mu <- 8
sigma <- 2
y <- rnorm(n=5,mean=mu,sd=sigma)
#Estimates
(ybar <- mean(y))
(sd <- sd(y))

#Calculate the density on a grid in the three cases
z.grid <- seq(mu-3*sigma,mu+3*sigma,length=1000)
d.true      <- dnorm(z.grid, mean=mu, sd=sigma)
d.plugin    <- dnorm(z.grid, mean=ybar, sd=sd)
d.sampleadj <- dTF(z.grid, mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1)

#Compute two-sided 95% PI
alpha <- 0.05
pi.true <- mu + qnorm(1-alpha/2) * c(-1,1) * sigma
pi.plugin <- ybar + qnorm(1-alpha/2) * c(-1,1) * sd
pi.sampleadj <- ybar + qt(1-alpha/2,df=n-1) * c(-1,1) * sd * sqrt(1+1/n)
mu.ci <- ybar + c(-1,1) * qnorm(1-alpha) * sd / sqrt(n)

tab <- rbind(pi.true, pi.plugin, pi.sampleadj, mu.ci)
print(tab,digits=2)

#How much more length?
(diff(pi.sampleadj)/diff(pi.plugin)-1)*100

#Compare the last one
pi <- qTF(c(alpha/2,1-alpha/2), mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1)

# Check coverage
diff(pTF(pi,  mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1))

#Coverage of different interval types
p_piplugin <- diff(pTF(pi.plugin, mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1))
p_muci <- diff(pTF(mu.ci, mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1))

## ----PLOTPREDINT--------------------------------------------------------------
layout(c(1,2),heights=c(5,1.5))
par(mar = c(2,5,1,1))

matplot(z.grid, cbind(d.plugin,d.sampleadj) ,type="l",xlab="z",ylab="density",col=c("royalblue","magenta"),lty=1,lwd=2)
rug(y,lwd=2)
legend(x="topleft",c("normal (plug-in)","non-standard t"),lwd=2,lty=1,col=c("royalblue","magenta"))

par(mar = c(2,5,1,1))

plot.pi <- function(pi,where=0,col="black") {
  lines( c(pi[1],pi[2]),rep(where,2),lty=1,lwd=2,col=col)
  lines( c(pi[1],pi[1]),where + c(-1,1)*0.25,lty=1,lwd=2,col=col)
  lines( c(pi[2],pi[2]),where + c(-1,1)*0.25,lty=1,lwd=2,col=col)
}

plot( c(pi.true[1],pi.true[2]),c(0,0),axes=FALSE,xlim=range(z.grid),xlab="",ylab="",ylim=c(0.75,2.25),col="orange",lty=1,lwd=2,type="n")

#plot.pi(pi.true,where=0,col="orange")
plot.pi(pi.plugin,where=1,col="royalblue")
plot.pi(pi.sampleadj,where=2,col="magenta")
#axis(1)
title(ylab="95% Pred.\n intervals",line=2)


## -----------------------------------------------------------------------------
#Setup monitoring and surveillance periods for the euromomo data
phase2 <- which(epoch(momo) >= "2007-10-01")
phase1 <- which(year(momo) == 2002 & epochInYear(momo) == 40):(phase2[1]-1)

######################################################################
# Illustrate how the Farrington algorithm works
######################################################################
par(mar=c(2, 4, 4, 2) + 0.1)
#Setup Farrington control object
alpha <- 0.005
cntrlFar <- list(range=phase2,alpha=alpha,b=5,w=4)#,powertrans="2/3")
#Adopt control object so exactly one computation is shown and
#setting the argument "plot" a picture is shown
one.t <- cntrlFar ; one.t$range <- min(one.t$range) ; one.t$plot <- TRUE
#Perform surveillance for one time point and save the graphics
onet <- farrington(momo[,"[75,85)"],control=one.t)


## ----echo=FALSE,results='hide', fig.keep="none"-------------------------------
#Load data and convert to an S4 object
data("hepatitisA")
hepatitisA <- disProg2sts(hepatitisA)

#Define parameters to use for Farrington algo
b <- 3
w <- 5
t0 <- 190

#Index of the reference values
t.ref <- 190 - (rep(1:b,each=2*w+1)*52 +  seq(-w,w,by=1))
#Make a data frame with the reference values
refvals <- data.frame(t=t.ref, y=observed(hepatitisA)[t.ref,])

#Show result
with(refvals,plot(t, y))

#Fit quasi-Poisson GLM
m <- glm(y ~ t, data=refvals, family=quasipoisson)
summary(m)

# Use the \texttt{predict} function to compute
#  $\hat{\mu}_{t_0}$ and $\sqrt{\Var(\hat{\mu}_{t_0})}$.
m.pred <- predict(m, newdata=data.frame(t=t0), se=TRUE,type="response")
phi <- summary(m)$dispersion

#Check computations
eta.hat <- coef(m) %*% c(1,t0)
mu.hat <- exp(eta.hat)
c(predict=m.pred$fit, manual=mu.hat)

var.log.mu.hat <- vcov(m)[1,1] + t0^2*vcov(m)[2,2] + 2*t0*vcov(m)[1,2]
c(predict=predict(m, newdata=data.frame(t=t0), se=TRUE,type="link")$se.fit,manual=sqrt(var.log.mu.hat))

var.mu.hat <- exp(eta.hat)^2 * var.log.mu.hat
c(predict.se=m.pred$se.fit,manual.se=sqrt(var.mu.hat))

#Aside: We actually know the true sd of the transformation,
#because result is a logN distribution
#sqrt(exp(var.log.mu.hat)-1)*exp(2*eta.hat + var.log.mu.hat)


#Compute the upper limit of a 99\% prediction interval for $y_{t_0}$
alpha <- 0.05
(pi <- m.pred$fit + c(-1,1)*qnorm(1-alpha/2)*sqrt( phi*m.pred$fit + m.pred$se.fit^2))

#Same result with surveillance package
control <- list(range=190,b=b,w=w,powertrans="none",alpha=alpha,trend=TRUE,reweight=FALSE,fitFun="algo.farrington.fitGLM")
surv <- farrington(hepatitisA, control)

upperbound(surv)


## -----------------------------------------------------------------------------
#Run the farrington algo and the improved version
cntrlFar$limit54 <- c(0,4)
cntrlFar2 <- modifyList(cntrlFar,
                        list(noPeriods=10,
                             populationOffset=FALSE,
                             fitFun="algo.farrington.fitGLM.flexible",
                             weightsThreshold=2.58,
                             pastWeeksNotIncluded=26,
                             pThresholdTrend=1,
                             thresholdMethod="nbPlugin"))

#Old call to farrington -- note that the prediction interval here is two sided
#and hence alpha has to be modified to be comparable
s.far <- farrington(momo[,"[75,85)"],control=modifyList(cntrlFar,list(alpha=alpha*2)))


#Show the results of the surveillance
plot(s.far,dx.upperbound=0,legend.opts=NULL,xlab="time (weeks)",alarm.symbol=list(pch=24,col=1,cex=1),col=c("darkgray",NA,1),lwd=c(1,1,2),lty=c(1,1,1),main="",ylim=c(-35,max(upperbound(s.far))))


#Add legend
legend(x=1,y=130,c("Alarm","NNBA"),pch=c(24,NA),lty=c(NA,1),horiz=TRUE,bg="white",lwd=c(2,3),col=c(1,1))


## ----eval=FALSE,echo=FALSE----------------------------------------------------
## alarmDates <- epoch(s.far[alarms(s.far) == 1,])
## formatDate(alarmDates,"%V")


## ----results="hide", message=FALSE--------------------------------------------
######################################################################
## R code generating graphs and numbers for the multivariate and scan
## statistics sections of the book chapter "Prospective Detection of
## Outbreaks" by B. Allévius and M. Höhle in the Handbook of
## Infectious Disease Data Analysis.
##
## Author: Benjamin Allévius <http://www.su.se/english/profiles/bekj9674-1.194276>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 2017-11-10
######################################################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(scales))

# Munge data ===================

# Menigococcal data from surveillance
data(imdepi)

# Aggregate across event types, age and gender
meningo_cases <- imdepi$events@data %>% as_tibble %>% mutate(count = 1L)

# Grab coordinates
meningo_coords <- coordinates(imdepi$events)
meningo_cases %<>%
  mutate(x_coord = meningo_coords[, 1], y_coord = meningo_coords[, 2])

# Add dates and remove unneeded columns
meningo_cases %<>%
  mutate(day = ceiling(time),
         date = as.Date("2001-12-31") + lubridate::days(day), # Exact start date is unknown
         year = as.integer(lubridate::year(date)),
         month = as.integer(month(date))) %>%
  dplyr::select(-eps.t, -eps.s, -.obsInfLength, -.sources, -.bdist, -.influenceRegion)

# Munge district-level data
district_data <- imdepi$stgrid %>% as_tibble %>%
  mutate(population = as.integer(popdensity * area))

tile_location <- tibble(tile = unique(district_data$tile),
                        location = 1:length(unique(district_data$tile)))
district_data <- left_join(district_data, tile_location, by = "tile")

# Get the total population (population is constant across time)
total_pop <- sum((district_data %>% filter(BLOCK == 1))$population)

# Monthly counts with covariates
meningo_monthly <- meningo_cases %>% dplyr::select(tile, BLOCK) %>%
  mutate(count = 1L) %>%
  group_by(tile, BLOCK) %>%
  summarize(count = n()) %>%
  ungroup %>%
  right_join(district_data %>% dplyr::select(tile, location, BLOCK, area, popdensity, population),
             by = c("tile", "BLOCK")) %>%
  mutate(count = ifelse(is.na(count), 0L, count)) %>%
  mutate(total_pop = total_pop) %>%
  rename(time = BLOCK) %>%
  arrange(location, time) %>%
  mutate(year = 2002 + floor((time - 1) / 12),
         month = ifelse(time %% 12 == 0, 12, time %% 12),
         date = as.Date(paste(year, month, "01", sep = "-")))



## ----message=FALSE------------------------------------------------------------
##Aggregate over space
meningo_ts <- meningo_monthly %>% group_by(year, month) %>%
  summarise(count=sum(count), population=sum(population),
            time=min(time),date=min(date)) %>%
  mutate(quarter=factor(paste0("Q",quarter(date)))) %>%
  ungroup


##Make a plot over time
ts <- ggplot( meningo_ts, aes(x=date, y=count)) +
  geom_crossbar(aes(ymin=0,ymax=count),fatten=0,fill="lightgray") +
  xlab("Time (month)") + ylab("No. cases") + ylim(c(0,NA)) +
  scale_x_date(date_breaks = "6 month",
               date_minor_breaks = "1 month",
               labels = scales::date_format("%b-%Y")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1))

##Look at it
ts


## ----warning=FALSE------------------------------------------------------------
# Scan statistics ==============================================================
suppressPackageStartupMessages(library(scanstatistics))

# Shapefile for the districts of Germany
load(system.file("shapes", "districtsD.RData", package = "surveillance"))

## ----compute_satscan, cache=TRUE----------------------------------------------
# Extract dates
dates <- (meningo_monthly %>% filter(tile == "01001") %>% dplyr::select(date))$date

# Parameters for the scan statistic
zones <- coordinates(districtsD) %>% coords_to_knn(k = 15) %>% knn_zones
scan_length <- 6 # Scanning window covers last 6 months

# Define surveillance period
t2_start <- 12 * 2 + 1
t2_end <- 12 * 4
t2_length <- t2_end - t2_start + 1

# Parameters for the surveillance period
scan_start <- t2_start
scan_end <- t2_end
scan_mc <- 999
scan_alpha <- 1 / (12 * 5)

# Store replicate scan statistics
replicates <- rep(NA, (scan_end - scan_start + 1) * scan_mc)

# Kulldorff's scan statistic
scan_df <- tibble(date = dates[scan_start:scan_end],
                  score = NA,
                  crit = NA,
                  pval = NA,
                  zone = NA,
                  duration = NA,
                  relrisk_in = NA,
                  relrisk_out = NA)

relrisk_support <- seq(1, 15, by = 0.1)
prev_relrisk_prob <- rep(1, length(relrisk_support))

# Run the scan statistics
idx <- 1
for (i in scan_start:scan_end) {
  time_window <- seq(max(1, i - scan_length + 1), i, by = 1)
  obs_counts <- meningo_monthly %>% filter(time %in% time_window)

  # Kulldorff's scan statistic
  scan <- scan_pb_poisson(obs_counts, zones, n_mcsim = scan_mc)
  repl_idx <- ((idx-1) * scan_mc + 1):(idx * scan_mc)
  replicates[repl_idx] <- scan$replicates$score

  scan_df$score[idx] <- scan$MLC$score
  scan_df$crit[idx] <- quantile(replicates[1:tail(repl_idx, 1)],
                                1 - scan_alpha,
                                type = 8)
  scan_df$pval[idx] <- scanstatistics:::mc_pvalue(scan$MLC$score, replicates[1:tail(repl_idx, 1)])
  scan_df$zone[idx] <- scan$MLC$zone_number
  scan_df$duration[idx] <- scan$MLC$duration
  scan_df$relrisk_in[idx] <- scan$MLC$relrisk_in
  scan_df$relrisk_out[idx] <- scan$MLC$relrisk_out

  print(paste0("idx = ", idx))
  idx <- idx + 1
}


## ----scan_score_plot, fig.align="center", warning=FALSE-----------------------
# Plot the score of the MLC over time
scan_score_plot <- ggplot(gather(scan_df,
                                 key = "type", value = "value",
                                 score, crit)) +
  geom_line(aes(x = date, y = value, color = type, linetype = type)) +
  scale_x_date(date_breaks = "6 month",
               date_minor_breaks = "1 month",
               labels = date_format("%b-%Y")) +
  scale_color_manual(name  = "",
                     breaks=c("score", "crit"),
                     labels=c("Score", "Critical value"),
                     values = c("gray47", "black")) +
  scale_linetype_manual(name  = "",
                        breaks=c("score", "crit"),
                        labels=c("Score", "Critical value"),
                        values = c("dashed", "solid")) +
  xlab("Date") + ylab(expression(lambda[W]^"*")) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank())
scan_score_plot


## ----scan_mlc_map, fig.align="center", warning=FALSE, message=FALSE-----------
# Calculate the overlap between clusters
zone_olap <- function(z1, z2) {
  length(base::intersect(z1, z2)) / length(base::union(z1, z2))
}

vec_zone_olap <- Vectorize(zone_olap, c("z1", "z2"))

zone_overlap <- rep(NA, nrow(scan_df) - 1)
for (i in 2:nrow(scan_df)) {
  zone_overlap[i-1] <- zone_olap(zones[[scan_df$zone[i - 1]]],
                                 zones[[scan_df$zone[i]]])
}

# Extract the MLC
MLC_zone <- zones[[scan_df$zone[which.max(scan_df$score)]]]

# Incidences per 100,000 people
meningo_incidence <- meningo_monthly %>%
  group_by(tile, location) %>%
  summarise(count = sum(count),
            population = population[1]) %>%
  ungroup %>%
  mutate(incidence = count * 100000 / population)

# Label the cluster
scan_clust <- districtsD@data %>%
  as_tibble %>%
  mutate(tile = as.factor(KEY),
         id = KEY,
         popdensity = POPULATION / AREA) %>%
  left_join(tile_location, by = "tile") %>%
  mutate(MLC = ifelse(location %in% MLC_zone, "Yes", "No")) %>%
  left_join(meningo_incidence, by = "tile")


# Make map data plotable
district_map <- fortify(districtsD) %>%
  as.tbl %>%
  left_join(scan_clust, by = "id") %>%
  mutate(state = substr(KEY, 1, 2))

cluster_state <- filter(district_map, MLC == "Yes" & order == 1)$state

district_map %<>%
  mutate(MLC_in_state = (state == cluster_state),
         MLC = factor(MLC, levels = c("Yes", "No")))


# Time series plot of observed scan statistic
scan_mlc_map <- ggplot(district_map %>% filter(MLC_in_state)) +
  theme_minimal() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = MLC),
               color = "black") +
  labs(x = "", y = "") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("gray37", "white")) +
  theme(legend.position = "none") +
  # theme(legend.position=c(0.9, 0.35)) +
  coord_equal()

library(ggspatial)
scan_mlc_map +
    annotation_north_arrow(which_north = "grid", location="br")

