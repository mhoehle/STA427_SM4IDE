## ----setup, include=FALSE, cache=FALSE----------------------------------------
#Load packages and initialize options
opts_chunk$set(fig.path='./knitr-figs/lecture03-',echo=FALSE, cache=TRUE,results='hide', warning=FALSE, fig.width=8,fig.height=4.5,size='scriptsize',error=FALSE,tidy=FALSE)
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
###  Rnw File for lecture 03
###
### History:
###  -- 2020-03-08 modified as part of STA427
###  -- 2020-06-22 file created as part of lecture 05 of SmitSjuk
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(lubridate)
library(RColorBrewer)


## -----------------------------------------------------------------------------
# Define time varying effective reproduction number
Ret <- function(date) ifelse(date <= as.Date("2020-03-15"), 2.5, 0.95)
# Define generation time to use
GT_pmf <- structure( c(0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1), names=0:7)
GT_obj <- R0::generation.time("empirical", val=GT_pmf)


## ----results="all"------------------------------------------------------------
GT_pmf


## ----simoutbreak, warning=FALSE-----------------------------------------------
####################################################################
#' Simulate time series with this generation time/serial interval and R_e(t)
#' @param n Number of time periods of the outbreak
#' @param Ret A function Re(t) which given t returns the current effective reproduc number
#' @param GT_obj
#' @param Initial number of cases for t=1
#' @return A names vector with the (expected) number of new symptom onsets per day. The outbreak is aligned such that $t=1$ corresponds to 2020-02-15.
####################################################################

routbreak <- function(n=100, Ret, GT_obj, initial_cases = 1) {
  # Set up time series of incident cases
  y <- rep(0, n + length(GT_pmf))
  y[seq_len(length(initial_cases))] <- initial_cases
  # Outbreak starts on 2020-02-15
  dates <- as.Date("2020-02-15") + 0:(n-1)
  # Extract serial interval PMF, ignore support at 0.
  GT_pmf <- GT_obj$GT[-1]

  # Loop over all time points
  for (i in 1:n) {
    date <- dates[i]
    y[i + 1:length(GT_pmf)] <- y[i] * (Ret(date) * GT_pmf) + y[i + 1:length(GT_pmf)]
  }

  # Data frame with the result. Assume we start on 15th of Feb
  res <- data.frame(Date=dates, y=y[1:n])

  #Done
  return(res)
}


# Generate an outbreak (no stochasticity, just the difference equation)
out <- routbreak(n=60, Ret=Ret, GT_obj=GT_obj)
out <- out %>% mutate(ratio = y/lag(y))
# Data frame with the true values
Ret_true <- data.frame(Date=out$Date) %>% mutate(R_hat=Ret(Date), Method="Truth")

## ----plotsimoutbreak, warning=FALSE-------------------------------------------
p1 <- ggplot(out, aes(x=Date, y=y)) + geom_line() + ylab("No. cases")
p2 <- ggplot(out, aes(x=Date, y=ratio)) + geom_line() + ylab("Ratio")
gridExtra::grid.arrange(p1, p2, ncol=2)


## ----epiestim, echo=TRUE, warning=FALSE---------------------------------------
library(EpiEstim)
# Rename data.frame columns to names handled by the EpiEstim pkg.
out_epiestim <- out %>% rename(I = y, dates = Date) %>% select(dates, I)

# Estimate the instantaneous reproduction number
res <- EpiEstim::estimate_R(out_epiestim, method = "non_parametric_si",
                            config=make_config(si_distr=GT_obj$GT,
                                               t_start=2:nrow(out_epiestim),
                                               t_end=2:nrow(out_epiestim))
)

# Convert result to a data.frame
rt_irt_df <- data.frame(Date=res$dates[res$R$t_end],
                        R_hat=res$R$`Mean(R)`,
                        lower=res$R$`Quantile.0.025`,
                        upper=res$R$`Quantile.0.975`,
                        Method="R(t)")


## ----PLOTINSTANTANEOUSR0, warning=FALSE, message=FALSE------------------------
ggplot(rt_irt_df, aes(x=Date, y=R_hat, color=Method)) +
  geom_ribbon(aes(x=Date,  ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  geom_line(data=Ret_true, aes(x=Date, y=R_hat, color=Method), lwd=2) +
  geom_line(lwd=1) +
  coord_cartesian(ylim=c(0, 5)) +
  scale_color_brewer(type="qual",palette="Set1", name="Method:") +
  ylab(expression(R(t))) +
  theme_minimal() +
  theme(legend.position="bottom")


## ----message=FALSE, warning=FALSE, results="hide"-----------------------------
#Load newest data from https://www.covid19.admin.ch/en/overview - for documentation
#see https://www.covid19.admin.ch/api/data/documentation
file_name <- file.path("Data", str_c("Swiss-COVID19-Data-",Sys.Date(),  ".zip"))
if (!file.exists(file_name)) {
  file_url <- "https://www.covid19.admin.ch/api/data/20210308-nsrvnhng/downloads/sources-csv.zip"
  download.file(url=file_url,destfile= file_name, mode="wb")
  utils::unzip(zipfile=file_name, exdir="Data")
}
# Load CSV file containing the case data for all regions
covid19_reports <- read_csv(file = file.path("Data", "data", "COVID19Cases_geoRegion.csv")) %>%
  filter(geoRegion == "CH")
# Hospitalizations
covid19_hosp <- read_csv(file = file.path("Data", "data", "COVID19Hosp_geoRegion.csv")) %>%
  filter(geoRegion == "CH")
# Deaths
covid19_deaths <- read_csv(file = file.path("Data", "data", "COVID19Death_geoRegion.csv")) %>%
  filter(geoRegion == "CH")

# I'll use values from
# https://www.drugs.com/medical-answers/covid-19-symptoms-progress-death-3536264/
# to define the lags
ts_df <- left_join(covid19_reports %>% rename(reported_cases = entries) %>% select(datum, reported_cases),
                   covid19_hosp %>% rename(reported_hosp = entries) %>% select(datum, reported_hosp), by=c("datum")) %>%
    left_join(covid19_deaths %>% rename(reported_deaths= entries) %>% select(datum, reported_deaths), by=c("datum")) %>%
  rename(Date = datum) %>%
  mutate(ratio_hosp_case = reported_hosp / lag(reported_cases, n=12)) %>% #12
  mutate(ratio_death_hosp = reported_deaths / lag(reported_hosp, n=7)) %>% #7
  mutate(ratio_death_case = reported_deaths / lag(reported_cases, n=19)) #19


# Helper function to calculate a centered running mean
roll7 <- function(x) { zoo::rollmean(x, k=7, align="center", na.pad=TRUE) }

# Compute running means for various quantities, replace this code by mutata_at
ts_df <- ts_df %>%
    mutate(reported_cases7 = roll7(reported_cases),
           reported_hosp7 = roll7(reported_hosp),
           reported_deaths7 = roll7(reported_deaths),
           ratio_hosp_case7 = roll7(ratio_hosp_case),
           ratio_death_hosp7 = roll7(ratio_death_hosp),
           ratio_death_case7 = roll7(ratio_death_case)
           )

## ----results="asis"-----------------------------------------------------------
# Show some summaries
ts_df %>%
  select(reported_cases, reported_hosp, reported_deaths) %>%
  summarise_all(sum, na.rm=TRUE) %>%
  xtable::xtable(digits=0) %>%
  print(include.rownames=FALSE)


## ----TS, warning=FALSE--------------------------------------------------------
## Show resulting time series
p1 <- ggplot(ts_df, aes(x=Date, y=reported_cases)) +
  geom_col() +
  geom_line(aes(y=reported_cases7), lwd=1.2)
p2 <- ggplot(ts_df, aes(x=Date, y=reported_hosp)) +
  geom_col() +
  geom_line(aes(y=reported_hosp7), lwd=1.2)
p3 <-ggplot(ts_df, aes(x=Date, y=reported_deaths)) +
  geom_col() +
  geom_line(aes(y=reported_deaths7), lwd=1.2)

## Show the three plots
gridExtra::grid.arrange(p1, p2, p3, ncol=1)


## ----warning=FALSE------------------------------------------------------------
# Plot of
p1 <- ggplot(ts_df, aes(x=Date, y=ratio_hosp_case)) + geom_line() +
  geom_line(aes(y=ratio_hosp_case7), col="steelblue", lwd=1.2) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"), NA), ylim=c(0,0.3)) +
  ggtitle("Hosp(t) / Cases(t-12)")
p2 <- ggplot(ts_df, aes(x=Date, y=ratio_death_hosp)) + geom_line() +
  geom_line(aes(y=ratio_death_hosp7), col="steelblue", lwd=1.2) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"), NA), ylim=c(0,1)) +
  ggtitle("Death(t) / Hosp(t-7)")
p3 <- ggplot(ts_df, aes(x=Date, y=ratio_death_case)) + geom_line() +
  geom_line(aes(y=ratio_death_case7), col="steelblue", lwd=1.2) +
  coord_cartesian(xlim=c(as.Date("2020-09-01"), NA), ylim=c(0,0.25)) +
  ggtitle("Death(t) / Cases(t-19)")

gridExtra::grid.arrange(p1, p2, p3, ncol=3)


## ----warning=FALSE, message=FALSE---------------------------------------------
#########################################
## R(t) estimation
#########################################

library(EpiEstim)

# Rename data.frame columns to names handled by the EpiEstim pkg.
out_epiestim <- ts_df %>% mutate(Date= as.Date(Date)) %>%
  rename(I = reported_cases, dates = Date) %>% select(dates, I) %>%
  # Ensure all reports are there (6 days back, this is rather conservative)
  filter(dates <= Sys.Date()-6)

# From the Swedish FOHM report: Här har vi använt det skattade serieintervallet från Nishiura et al.
# (2020) med ett medelvärde på 4.8 dagar och en standaravvikelse på 2.3 dagar.
# Note: Weibull distribution in the Nishiura et al. paper vs. (shifted) gamma in EpiEstim
# Note sure if BAG uses an own estimate?

# Estimate the instantaneous R_e(t) - no smoothing
res1 <- EpiEstim::estimate_R(out_epiestim, method = "parametric_si",
                            config=make_config(mean_si = 4.8, std_si = 2.3,
                                               t_start=2:nrow(out_epiestim),
                                               t_end=2:nrow(out_epiestim)))
# Use smoothing over 7 days (recommended!)
res7 <- EpiEstim::estimate_R(out_epiestim, method = "parametric_si",
                            config=make_config(mean_si = 4.8, std_si = 2.3,
                                               t_start=2:(nrow(out_epiestim)-6),
                                               t_end=8:nrow(out_epiestim)))

# Convert result to a data.frame
epiestim2df <- function(res) {
  data.frame(Date=res$dates[res$R$t_end], R_hat=res$R$`Mean(R)`, lower=res$R$`Quantile.0.025(R)`, upper=res$R$`Quantile.0.975`)
}

rt_irt7_df <- epiestim2df(res7) %>% mutate(method = "smooth7")
rt_irt1_df <- epiestim2df(res1) %>% mutate(method = "no smoothing")
rt_df <- rbind(rt_irt7_df, rt_irt1_df)

#Use R0 pkg
# Define generation time to use
GT_obj <- R0::generation.time("lognormal", val=c( 4.8, 2.3))

# Use Wallinga-Teunis estimator
ts <- out_epiestim %>% mutate(Date=row_number(), y=I)
est <- R0::est.R0.TD(ts$y, GT_obj, begin=1L, end=nrow(ts), nsim=100)
wt <- out_epiestim %>% rename(Date = dates) %>%
  mutate(R_hat = est$R, lower=NA, upper=NA)

#Read BAG estimate from the files
rt_bag <- read_csv(file.path("Data","data","COVID19Re_geoRegion.csv")) %>%
  filter(geoRegion == "CH") %>%
  mutate(method ="smooth7", type="BAG",lower=NA, upper=NA) %>%
  rename(Date=date, R_hat=median_R_mean) %>%
  select(Date, R_hat, lower, upper,method, type)

## ----rtplots, warning=FALSE---------------------------------------------------
# The plot - doesn't make too much sense because testing is substantially
# increased in the last weeks
# Start at 2020-03-16 so imported casesd do not play a role anymore
ggplot(rt_df, aes(x=Date, y=R_hat)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="steelblue", alpha=0.2) +
  geom_line(aes(color="EpiEstim")) +
  geom_line(data=rt_bag, aes(color="BAG estimate"),lty=2) +
  #geom_line(data=wt, aes(color="Wallinga-Teunis"),lty=2) +
  xlab("Date (of report)") +
  #scale_color_manual(values=c("BAG estimate"="#7FC97F", "EpiEstim"="#000000", "Wallinga-Teunis"="#BEAED4"), name="") +
  scale_color_manual(values=c("BAG estimate"="#7FC97F", "EpiEstim"="#000000"), name="") +
  coord_cartesian(xlim=c(as.Date("2020-03-01"), NA), ylim=c(0, 2)) +
  ylab(expression(R(t)))  +
  geom_hline(yintercept=1, lty=2, col="salmon2") +
  facet_wrap(~ method, ncol=1)


## ----warning=FALSE, message=FALSE---------------------------------------------
covid19_tests <- read_csv(file = file.path("Data", "data", "COVID19Test_geoRegion_all.csv")) %>%
  filter(geoRegion == "CH") %>% rename(Date=datum)

covid19_tests_long <- covid19_tests %>%
  mutate(proportion_pos = entries_pos / entries) %>%
  select(Date, entries, proportion_pos) %>%
  mutate(entries7 = roll7(entries),
         proportion_pos7 = roll7(proportion_pos)) %>%
  pivot_longer(-Date, names_to = "type", values_to = "value") %>%
  mutate( type2 =  case_when( type == "entries"  ~  "Total no. test",
                              type == "entries7" ~  "Total no. test",
                             type == "proportion_pos" ~ "Proportion positive",
                             type == "proportion_pos7" ~ "Proportion positive"),
          smooth7 = str_detect(type, "7"))


ggplot(covid19_tests_long, aes(x=Date, y=value, color=smooth7)) +
  geom_line() +
  ylab("Proportion / No. tests") +
  scale_color_brewer(palette="Paired") +
  coord_cartesian(xlim=(c(ymd("2020-06-01"),NA))) + #collection doesnt start before july
  facet_wrap(~ type2, ncol=1, scale="free_y")

