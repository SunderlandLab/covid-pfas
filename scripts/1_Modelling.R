################################################################################
# Author: Jahred M. Liddie (reviewed by MK)
# Purpose: Primary models for associations between PFAS DW contamination 
  # and COVID-19 cumulative mortality
# Date began: 2/22/23
################################################################################
library(tidyverse)
library(MASS)
library(lme4)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load 'statewide' and UCMR processed data
dat <- read.csv("../data/processed/covid_CWS.csv")

dat.ucmr <- read.csv("../data/processed/covid_UCMR.csv")

# total deaths by end of study period in each dataset
sum(round(dat$deaths, -2))
sum(round(dat.ucmr$deaths, -2))

################################################################################
# Modelling with state sampling data
################################################################################
# unadjusted model
b1 <- glm(deaths ~ Any_bin, 
          offset = log(tot_pop), 
          family = "quasipoisson", data = dat)

c1 <- glm(deaths ~ Any_reg_bin, 
          offset = log(tot_pop), 
          family = "quasipoisson", data = dat)

# adding all other potential confounders and predictors
n1a <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                offset(log(tot_pop)), 
                random = ~ 1 | state,
                family = "quasipoisson", data = dat)

p1 <- glmmPQL(deaths ~ Any_reg_bin + scale(days_since_fc) + scale(bed.1000) + 
                 scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                 scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                 scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                 offset(log(tot_pop)), 
               random = ~ 1 | state,
               family = "quasipoisson", data = dat)

################################################################################
# Modelling with UCMR 3 data
################################################################################
# unadjusted model
d1 <- glm(deaths ~ Any5, 
          offset = log(tot_pop), 
          family = "quasipoisson", data = dat.ucmr)

# now add all other covariates
m1 <- glmmPQL(deaths ~ Any5 + days_since_fc + scale(bed.1000) + 
                scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                scale(log(median_income)) + offset(log(tot_pop)), 
                random = ~ 1 | state,
                family = "quasipoisson", data = dat.ucmr)

################################################################################
# now extract effect estimates from final models:
options(digits = 3)

extract_CIR.f <- function(model = NULL, exposure.desc = NULL) {
  
  model.summary <- summary(model)
  model.class <- class(model)[1]
  
  if (model.class == "glm") {
    
    model.results <- tibble(
      logCIR = model.summary$coefficients[2],
      logSE = model.summary$coefficients[2,2],
      nobs = model.summary$df.residual + 2, # nobs calculated for this unadj glm
      exposure.desc = exposure.desc,
      CIR = exp(logCIR),
      CIR.LCI = exp(logCIR - 1.96*logSE),
      CIR.UCI = exp(logCIR + 1.96*logSE),
      full.estimate = paste(round(CIR,2), " [95% CI: ", round(CIR.LCI,2), ", ", round(CIR.UCI,2), "]", sep = "")
    )
    
  }


  else{
    
    model.results <- tibble(
      logCIR = model.summary$coefficients$fixed[2],
      logSE = sqrt( model.summary$varFix[2,2] ), # standard error
      nobs = nrow(model.summary$data),
      exposure.desc = exposure.desc,
      CIR = exp(logCIR),
      CIR.LCI = exp(logCIR - 1.96*logSE),
      CIR.UCI = exp(logCIR + 1.96*logSE),
      full.estimate = paste(round(CIR,2), " [95% CI: ", round(CIR.LCI,2), ", ", round(CIR.UCI,2), "]", sep = "")
    )

  }
  
  return(model.results)
}

# final models are: n1a, p1, o1, and m1
n1a_results <- extract_CIR.f(n1a, exposure.desc = "Any_bin (Statewide)")
p1_results <- extract_CIR.f(p1, exposure.desc = "Any_reg (Statewide)")
m1_results <- extract_CIR.f(m1, exposure.desc = "Any_bin (UCMR 3)")

unadj1 <- extract_CIR.f(b1, exposure.desc = "Any_bin (Statewide; unadjusted)")
unadj2 <- extract_CIR.f(c1, exposure.desc = "Any_reg (Statewide; unadjusted)")
unadj3 <- extract_CIR.f(d1, exposure.desc = "Any_bin (UCMR 3; unadjusted)")

all.results <- rbind(
  n1a_results, p1_results, m1_results, unadj1, unadj2, unadj3
)

if (FALSE) {
  write.csv(all.results, "../results/main_regression_results.csv")
}
