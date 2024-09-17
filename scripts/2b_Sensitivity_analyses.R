################################################################################
# Author: Jahred Liddie (reviewed by MK)
# Purpose: Sensitivity analyses for associations between PFAS DW contamination 
  # and COVID-19 cumulative mortality
# Date began: 2/22/23
################################################################################
library(tidyverse)
library(MASS)
library(lme4)
library(glmmTMB)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load 'statewide' and UCMR processed data; exclude counties without all PFAS data
dat <- read.csv("../data/processed/covid_CWS.csv")

dat.ucmr <- read.csv("../data/processed/covid_UCMR.csv")

################################################################################
# data-focused sensitivity analyses
# 1. excluding CA, PA, KY, NY, SC, UT, and MD (states with any sort of targeted sampling)
states_excl <- c("CA", "PA", "KY", "UT", "MD")

n1_exclstates <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                           scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                           scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                           scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                           scale(log(median_income)) + offset(log(tot_pop)),
                         random = ~ 1 | state,
                         family = "quasipoisson", data = subset(dat, ! (state %in% states_excl)))

# 2. excluding LA, Cook County IL, and New York metropolitan area
outlying_areas <- c("ca_los angeles", "il_cook", "ny_new york")
dat_nooutliers <- dat %>% filter(!state_county %in% outlying_areas)

n1_nooutliers <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                           scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                           scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                           scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                           scale(log(median_income)) + offset(log(tot_pop)),
                         random = ~ 1 | state,
                         family = "quasipoisson", data = dat_nooutliers)

dat.ucmr_nooutliers <- dat.ucmr %>% filter(!state_county %in% outlying_areas)

m1_nooutliers <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                          scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                          scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                          scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                          scale(log(median_income)) + offset(log(tot_pop)), 
                        random = ~ 1 | state,
                        family = "quasipoisson", data = dat.ucmr_nooutliers)

################################################################################
# exposure-focused sensitivity analyess
# 1. fraction of systems with contamination
n1_frac <- glmmPQL(deaths ~ scale(frac_Any_detect) + scale(days_since_fc) + scale(bed.1000) + 
                       scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                       scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                       scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                       scale(log(median_income)) + offset(log(tot_pop)),
                     random = ~ 1 | state,
                     family = "quasipoisson", data = dat)

sd(dat$frac_Any_detect)

# 2. fraction of systems w/ contamination weighted by population served
n1_weight <- glmmPQL(deaths ~ scale(frac_Any_detect_w) + scale(days_since_fc) + scale(bed.1000) + 
                       scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                       scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                       scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                       scale(log(median_income)) + offset(log(tot_pop)),
                       random = ~ 1 | state,
                     family = "quasipoisson", data = dat)

sd(dat$frac_Any_detect_w)

# 3. including PFHpA (UCMR 3)
m1_PFHpA <- glmmPQL(deaths ~ Any6 + scale(days_since_fc) + scale(bed.1000) + 
                     scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                     scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                     scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                     scale(log(median_income)) + offset(log(tot_pop)), 
                   random = ~ 1 | state,
                   family = "quasipoisson", data = dat.ucmr)

################################################################################
# outcome-focused sensitivity analyses
# 1. 2/3rds of way to end-date
n1_threeq <- glmmPQL(deaths_threeq ~ Any_bin + scale(days_since_fc_threeq) + scale(bed.1000) + 
                           scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                           scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                           scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                           scale(log(median_income)) + offset(log(tot_pop)),
                         random = ~ 1 | state,
                         family = "quasipoisson", data = dat)

m1_threeq <- glmmPQL(deaths_threeq ~ Any5 + scale(days_since_fc_threeq) + scale(bed.1000) + 
                       scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                       scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                       scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                       scale(log(median_income)) + offset(log(tot_pop)), 
                     random = ~ 1 | state,
                     family = "quasipoisson", data = dat.ucmr)

# 2. instead using an offset for the number of ppl over age 65
n1_65offset <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                         scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                         scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                         scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                         scale(log(median_income)) + offset(log(pop65plus)),
                       random = ~ 1 | state,
                       family = "quasipoisson", data = dat)

m1_65offset <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                         scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                         scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                         scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                         scale(log(median_income)) + offset(log(pop65plus)),
                       random = ~ 1 | state,
                       family = "quasipoisson", data = dat.ucmr)

################################################################################
# modelling-focused sensitivity analyses (or additional adjustments)
# 1. fixed effects
n1_FEs <- glm(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                     scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                     scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                     scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                     state,
                   offset = log(tot_pop), 
                   family = "quasipoisson", data = dat)

m1_FEs <- glm(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                state,
              offset = log(tot_pop), 
              family = "quasipoisson", data = dat.ucmr)

# 2. Mixed model using county-level intercepts (OLRE model)
n1_county <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                         scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                         scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                         scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                         scale(log(median_income)) + offset(log(pop65plus)),
                       random = ~ 1 | state_county,
                       family = "quasipoisson", data = dat)

m1_county <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                       scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                       scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                       scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                       scale(log(median_income)) + offset(log(tot_pop)), 
                     random = ~ 1 | state_county,
                     family = "quasipoisson", data = dat.ucmr)
  
# 3. negative binomial
n1_negbin <- glmmTMB::glmmTMB(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                                scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                                scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                                scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                                scale(log(median_income)) + offset(log(tot_pop)) + (1 | state),
                              family = nbinom1(),
                              data = dat)

m1_negbin <- glmmTMB::glmmTMB(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                              scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                              scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                              scale(log(med_OOHH_value)) + scale(percunderHSed) + 
                              scale(log(median_income)) + offset(log(tot_pop)) + (1 | state),
                              family = nbinom1(),
                            data = dat.ucmr)
  
# 4. adding a proxy variable for susceptibility prior to the pandemic: 
  # county-level all-cause mortality rate in 2015
n1_prior <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                      scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                      scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                      scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                      scale(all_cause_rate_2015) + offset(log(tot_pop)),
                    random = ~ 1 | state,
                    family = "quasipoisson", data = dat)

m1_prior <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                      scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                      scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                      scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                      scale(all_cause_rate_2015) + offset(log(tot_pop)),
                    random = ~ 1 | state,
                    family = "quasipoisson", data = dat.ucmr)

# 4. control for number of tests
n1_tests <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                      scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                      scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                      scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                      scale(totalTestResults) + offset(log(tot_pop)),
                    random = ~ 1 | state,
                    family = "quasipoisson", data = dat)

m1_tests <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                      scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                      scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                      scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                      scale(totalTestResults) + offset(log(tot_pop)),
                    random = ~ 1 | state,
                    family = "quasipoisson", data = dat.ucmr)

# 5. control a function of lat + long
n1_latlong <- mgcv::gamm(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                          scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                          scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                          scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                          s(lat, bs = "cs") + s(lon, bs = "cs") + offset(log(tot_pop)),
                         random = list(state=~1),
                        family = "quasipoisson", data = dat)

m1_latlong <- mgcv::gamm(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                          scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                          scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                          scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                          s(lat, bs = "cs") + s(lon, bs = "cs") + offset(log(tot_pop)),
                        random = list(state=~1),
                        family = "quasipoisson", data = dat.ucmr)

# 5. control for adult smoking
n1_smoking <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) +
                        scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                        scale(smoking) + offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat)

m1_smoking <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) +
                      scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                      scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                      scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                      scale(smoking) + offset(log(tot_pop)),
                    random = ~ 1 | state,
                    family = "quasipoisson", data = dat.ucmr)

# 5. control for pm2.5
n1_pm25 <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) +
                        scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                        scale(pm25) + offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat)

m1_pm25 <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) +
                        scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                        scale(pm25) + offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat.ucmr)

# 5. control for MCL violations
n1_MCL <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) +
                     scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                     scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                     scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                     scale(MCL_hbv_nc) + offset(log(tot_pop)),
                   random = ~ 1 | state,
                   family = "quasipoisson", data = dat)

m1_MCL <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) +
                     scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                     scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                     scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                     scale(MCL_hbv_nc) + offset(log(tot_pop)),
                   random = ~ 1 | state,
                   family = "quasipoisson", data = dat.ucmr)

# 5. drop median OOHH value
  # drop1 <- glmmPQL(deaths ~ scale(log(med_OOHH_value)) + offset(log(tot_pop)),
  #                  random = ~ 1 | state,
  #                  family = "quasipoisson",
  #                  data = dat)
  # 
  # drop2 <- glmmPQL(deaths ~ scale(log(median_income)) + offset(log(tot_pop)),
  #                  random = ~ 1 | state,
  #                  family = "quasipoisson",
  #                  data = dat)
  # # coefficient for median_income is larger, so dropping OOHH value:

n1_drop1 <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) +
                    scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                    scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                    scale(percunderHSed) + scale(log(median_income)) +
                    offset(log(tot_pop)),
                  random = ~ 1 | state,
                  family = "quasipoisson", data = dat)

m1_drop1 <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) +
                    scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                    scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                    scale(percunderHSed) + scale(log(median_income)) +
                    offset(log(tot_pop)),
                  random = ~ 1 | state,
                  family = "quasipoisson", data = dat.ucmr)

# 6. drop hospital beds
  # drop3 <- glmmPQL(deaths ~ scale(pop.dens.sqmi) + offset(log(tot_pop)),
  #                  random = ~ 1 | state,
  #                  family = "quasipoisson",
  #                  data = dat)
  # 
  # drop4 <- glmmPQL(deaths ~ scale(bed.1000) + offset(log(tot_pop)),
  #                  random = ~ 1 | state,
  #                  family = "quasipoisson",
  #                  data = dat)
  # # coefficient for population density is larger, so dropping hospital beds:

n1_drop2 <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + 
                    scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                    scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                    scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                    offset(log(tot_pop)),
                  random = ~ 1 | state,
                  family = "quasipoisson", data = dat)

m1_drop2 <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + 
                    scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                    scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                    scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) +
                    offset(log(tot_pop)),
                  random = ~ 1 | state,
                  family = "quasipoisson", data = dat.ucmr)

################################################################################
# confounding-focused sensitivity analyses (or additional adjustments)
# 1. vary population density adjustments (quartiles and log)
dat$pop.dens.sqmi_quart <- cut(dat$pop.dens.sqmi, 
                               quantile(dat$pop.dens.sqmi, seq(0, 1, by = 0.25)), 
                               include.lowest = TRUE)

dat.ucmr$pop.dens.sqmi_quart <- cut(dat.ucmr$pop.dens.sqmi, 
                                   quantile(dat.ucmr$pop.dens.sqmi, seq(0, 1, by = 0.25)), 
                                   include.lowest = TRUE)

n1_popdensquart <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                        pop.dens.sqmi_quart + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                        offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat)

n1_logpopdens <- glmmPQL(deaths ~ Any_bin + scale(days_since_fc) + scale(bed.1000) + 
                             log(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                             scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                             scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                             offset(log(tot_pop)),
                           random = ~ 1 | state,
                           family = "quasipoisson", data = dat)

m1_popdensquart <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                        pop.dens.sqmi_quart + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                        offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat.ucmr)

m1_logpopdens <- glmmPQL(deaths ~ Any5 + scale(days_since_fc) + scale(bed.1000) + 
                        log(pop.dens.sqmi)  + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                        offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat.ucmr)

# 2. vary control of days since first case
dat$days_since_fc_quart <- cut(dat$days_since_fc, 
                               quantile(dat$days_since_fc, seq(0, 1, by = 0.25)), 
                               include.lowest = TRUE)

dat.ucmr$days_since_fc_quart <- cut(dat.ucmr$days_since_fc, 
                                    quantile(dat.ucmr$days_since_fc, seq(0, 1, by = 0.25)), 
                                    include.lowest = TRUE)

n1_dsfc <- glmmPQL(deaths ~ Any_bin + days_since_fc_quart + scale(bed.1000) + 
                     pop.dens.sqmi_quart + scale(percHisp) + scale(percWhite) +
                     scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                     scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                     offset(log(tot_pop)),
                   random = ~ 1 | state,
                   family = "quasipoisson", data = dat)

m1_dsfc <- glmmPQL(deaths ~ Any5 + days_since_fc_quart + scale(bed.1000) + 
                     pop.dens.sqmi_quart + scale(percHisp) + scale(percWhite) +
                     scale(percBlack) + scale(perc65plus) + scale(perchomeowner) +
                     scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                     offset(log(tot_pop)),
                   random = ~ 1 | state,
                   family = "quasipoisson", data = dat.ucmr)

# 4. control for age categories differently 
dat$perc65plus_quart <- cut(dat$perc65plus, 
                               quantile(dat$perc65plus, seq(0, 1, by = 0.25)), 
                               include.lowest = TRUE)

dat.ucmr$perc65plus_quart <- cut(dat.ucmr$perc65plus, 
                                    quantile(dat.ucmr$perc65plus, seq(0, 1, by = 0.25)), 
                                    include.lowest = TRUE)

n1_p65_quart <- glmmPQL(deaths ~ Any_bin + days_since_fc + scale(bed.1000) + 
                         scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                         scale(percBlack) + perc65plus_quart + scale(perchomeowner) +
                         scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                         offset(log(tot_pop)),
                       random = ~ 1 | state,
                       family = "quasipoisson", data = dat)

m1_p65_quart <- glmmPQL(deaths ~ Any5 + days_since_fc + scale(bed.1000) + 
                          scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                          scale(percBlack) + perc65plus_quart + scale(perchomeowner) +
                          scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                          offset(log(tot_pop)),
                        random = ~ 1 | state,
                        family = "quasipoisson", data = dat.ucmr)

# 5. control for additional age categories as continuous vars
n1_moreage <- glmmPQL(deaths ~ Any_bin + days_since_fc + scale(bed.1000) + 
                          scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                          scale(percBlack) + scale(perc65plus) + scale(perc45_64) +
                          scale(perc15_44) + scale(perchomeowner) +
                          scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                          offset(log(tot_pop)),
                        random = ~ 1 | state,
                        family = "quasipoisson", data = dat)

m1_moreage <- glmmPQL(deaths ~ Any5 + days_since_fc + scale(bed.1000) + 
                        scale(pop.dens.sqmi) + scale(percHisp) + scale(percWhite) +
                        scale(percBlack) + scale(perc65plus) + scale(perc45_64) +
                        scale(perc15_44) + scale(perchomeowner) +
                        scale(log(med_OOHH_value)) + scale(percunderHSed) + scale(log(median_income)) + 
                        offset(log(tot_pop)),
                      random = ~ 1 | state,
                      family = "quasipoisson", data = dat.ucmr)

# ##############################################################################
# Combine all estimates into one df
# ##############################################################################
# now extract effect estimates from final models:
options(digits = 3)

extract_CIR.f <- function(model = NULL, type = "binary", exposure.desc = NULL) {
  
  model.type <- class( get(model) )[1]
  
  if (model.type == "gamm") {
    
    model.summary <- summary(get(model)$lme)
  
  }
  
  else{
    
  model.summary <- summary( get(model) ) # need to "get" model since it will be characters
  
  }
  
  model.results <- tibble(
      logCIR = ifelse(grepl("^glmmPQL$", model.type) | grepl("^gamm$", model.type), 
                      model.summary$coefficients$fixed[2], 
                      ifelse(grepl("^glmmTMB$", model.type), model.summary$coefficients$cond[2,1],
                             model.summary$coefficients[2])
                      ),
      
      logSE = ifelse(grepl("^glmmPQL$", model.type) | grepl("^gamm$", model.type),
                     sqrt( model.summary$varFix[2,2] ), 
                     ifelse(grepl("^glmmTMB$", model.type), model.summary$coefficients$cond[2,2],
                            model.summary$coefficients[2,2])
                     ),
      
      nobs = ifelse(grepl("^glmmTMB$", model.type), model.summary$nobs,
                           ifelse(grepl("^gamm$", model.type), nrow(get(model)$lme$data),
                                  nrow(get(model)$data)
                                  )
                           ),
      
      CIR = exp(logCIR),
      CIR.LCI = exp(logCIR - 1.96*logSE),
      CIR.UCI = exp(logCIR + 1.96*logSE),
      model = paste(model),
      full.estimate = paste(round(CIR,2), " [", round(CIR.LCI,2), ", ", round(CIR.UCI,2), "]", sep = ""),
      class = model.type
    )
  
  return(model.results)

}

all.objects <- ls()
statewide_models <- c(all.objects[grepl("^n1_", all.objects)])

ucmr_models <- c(all.objects[grepl("^m1_", all.objects)])

statewide_results <- map_dfr(statewide_models,
                                ~extract_CIR.f(.x, type = "binary", 
                                               exposure.desc = c("Statewide data")))

statewide_results <- statewide_results %>%
  mutate(category = case_when(grepl("outlier", model) | grepl("excl", model) ~ "Data",
                              grepl("frac", model) | grepl("weight", model) ~ "Exposure",
                              grepl("offset", model) | grepl("threeq", model) ~ "Outcome",
                              grepl("FE", model) | grepl("negbin", model) | grepl("county", model) ~ "Modelling",
                              grepl("latlong", model) | grepl("p65", model) | grepl("pm", model) | grepl("smok", model) |
                              grepl("popdens", model) | grepl("prior", model) | grepl("tests", model) |
                              grepl("dsfc", model) | grepl("MCL", model) | grepl("age", model) | grepl("drop", model) ~ "Covariates",
                              TRUE ~ NA),
         dataset = "Statewide sampling",
         sens_name = case_when(grepl("outlier", model) ~ "Excl. outliers",
                               grepl("excl", model) ~ "Excl. certain states",
                               grepl("frac", model) ~ "% systems detecting",
                               grepl("weight", model) ~ "Weighted % systems detecting",
                               grepl("offset", model) ~ "Offset: total pop. 65+",
                               grepl("threeq", model) ~ "Different end-date",
                               grepl("FE", model) ~ "State fixed-effects",
                               grepl("negbin", model) ~ "Negative binomial mixed model",
                               grepl("county", model) ~ "County random-effects",
                               grepl("latlong", model) ~ "Incl. lat + long",
                               grepl("p65", model) ~ "Quartiles of percent 65+", 
                               grepl("age", model) ~ "Add. age categories", 
                               grepl("pm", model) ~ "Incl. long-term PM2.5",
                               grepl("smok", model) ~ "Incl. % adult smokers",
                               grepl("popdensquart", model) ~ "Quartiles of pop. density",
                               grepl("logpopdens", model) ~ "Log pop. density",
                               grepl("prior", model) ~ "Incl. baseline mortality rate",
                               grepl("tests", model) ~ "Incl. total # of tests",
                               grepl("dsfc", model) ~ "Quartiles of days since first case",
                               grepl("MCL", model) ~ "Incl. prior MCL violations",
                               grepl("drop1", model) ~ "Excl. median OOHH value",
                               grepl("drop2", model) ~ "Excl. hospital beds",
                               TRUE ~ NA),
         )

ucmr_results <- map_dfr(ucmr_models,
                        ~extract_CIR.f(.x, type = "binary", 
                                       exposure.desc = c("UCMR 3")))

ucmr_results <- ucmr_results %>%
  mutate(category = case_when(grepl("outlier", model) | grepl("excl", model) ~ "Data",
                              grepl("PFHpA", model) ~ "Exposure",
                              grepl("offset", model) | grepl("threeq", model) ~ "Outcome",
                              grepl("FE", model) | grepl("negbin", model) | grepl("county", model) ~ "Modelling",
                              grepl("race", model) | grepl("latlong", model) | grepl("p65", model) |
                              grepl("popdens", model) | grepl("prior", model) | grepl("tests", model) |
                              grepl("dsfc", model) | grepl("smok", model) | grepl("pm", model) |
                              grepl("MCL", model) | grepl("age", model) | grepl("drop", model) ~ "Covariates",
                              TRUE ~ NA),
         dataset = "UCMR 3",
         sens_name = case_when(grepl("outlier", model) ~ "Excl. outliers",
                               grepl("excl", model) ~ "Excl. certain states",
                               grepl("PFHpA", model) ~ "Incl. PFHpA",
                               grepl("offset", model) ~ "Offset: total pop. 65+",
                               grepl("threeq", model) ~ "Different end-date",
                               grepl("FE", model) ~ "State fixed-effects",
                               grepl("county", model) ~ "County random-effects",
                               grepl("negbin", model) ~ "Negative binomial mixed model",
                               grepl("latlong", model) ~ "Incl. lat + long",
                               grepl("p65", model) ~ "Quartiles of percent 65+",
                               grepl("age", model) ~ "Add. age categories", 
                               grepl("pm", model) ~ "Incl. long-term PM2.5",
                               grepl("smok", model) ~ "Incl. % adult smokers",
                               grepl("popdensquart", model) ~ "Quartiles of pop. density",
                               grepl("logpopdens", model) ~ "Log pop. density",
                               grepl("prior", model) ~ "Incl. baseline mortality rate",
                               grepl("tests", model) ~ "Incl. total # of tests",
                               grepl("dsfc", model) ~ "Quartiles of days since first case",
                               grepl("MCL", model) ~ "Incl. prior MCL violations",
                               grepl("drop1", model) ~ "Excl. median OOHH value",
                               grepl("drop2", model) ~ "Excl. hospital beds",
                               TRUE ~ NA),
  )

all_sensitivity <- rbind(statewide_results, ucmr_results)

# then save here:
if (FALSE){
  write.csv(all_sensitivity, "../results/all_sensitivity.csv")
}
