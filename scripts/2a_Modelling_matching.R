################################################################################
# Author: Jahred Liddie (reviewed by MK)
# Purpose: Matching methods as sensitivity analysis/robustness check
# date created: 2/22/23
################################################################################
library(tidyverse)
library(MatchIt)
library(cobalt)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load statewide and UCMR 3 processed data
dat.statewide <- read.csv("../data/processed/covid_CWS.csv")

# selecting primary covariates
dat.statewide <- dat.statewide %>%
  mutate(census_region = as.factor(census_region),
         state = as.factor(state)) %>%
  dplyr::select(state_county, Any_bin, Any_reg_bin, census_region, state, percBlack, 
                percWhite, percHisp, percunderHSed, median_income, 
                med_OOHH_value, perchomeowner, perc65plus, pop.dens.sqmi, 
                days_since_fc, bed.1000, deaths, tot_pop)

dat.ucmr <- read.csv("../data/processed/covid_UCMR.csv")

dat.ucmr <- dat.ucmr %>%
  mutate(census_region = as.factor(census_region),
         state = as.factor(state)) %>%
  dplyr::select(state_county, state, census_region, 
                Any5, percBlack, percWhite, percHisp, percunderHSed,
                median_income, med_OOHH_value, perchomeowner, perc65plus, 
                pop.dens.sqmi, tot_pop, days_since_fc, bed.1000, deaths)

################################################################################
# function to extract results from MatchIt
extract_results.f <- function(matched.model = NULL, matched.dataset = NULL, 
                              full.dataset = NULL, match.type = NULL,
                              model.adjust = NULL) {

  estimate <- exp( coef(matched.model) )[2]
  
  LCL <- exp( confint( matched.model ) )[2,1]
  
  UCL <- exp( confint( matched.model ) )[2,2]
  
  print(LCL)
  print(UCL)
  
  # sample sizes
  matches <- n_distinct(matched.dataset$subclass)
  n.match <- nrow(matched.dataset)
  n.total <- nrow(full.dataset)
  
  results <- tibble(
    estimate = estimate,
    LCL = LCL,
    UCL = UCL,
    matches = matches,
    n.match = n.match,
    n.total = n.total,
    p.matched = n.match/n.total,
    match.type = match.type,
    model.adjust = model.adjust
  )
  
  return(results)
}

################################################################################
# Propensity score matching (nearest 1:1);
  # Notes: here, the caliper is set to 1, matching without replacement, 
  # and propensity scores are re-estimated after discarding
  # treated and untreated units
m_statewide_PSM <- matchit(
  Any_bin ~ bed.1000 + percHisp + pop.dens.sqmi +
    percBlack + perc65plus + perchomeowner +
    log(med_OOHH_value) + percunderHSed + percWhite +
    log(median_income) + days_since_fc,
  data = dat.statewide,
  method = "nearest",
  distance = "glm",
  exact = ~ census_region,
  caliper = 1,
  discard = "both",
  replace = FALSE,
  reestimate = TRUE)

match_statewide_PSM <- match.data(m_statewide_PSM, distance = "prop.score")
summary(m_statewide_PSM)

statewide_PSM <- glm(deaths ~ Any_bin, offset = log(tot_pop),
                          family = "quasipoisson", data = match_statewide_PSM)

r1 <- extract_results.f(matched.model = statewide_PSM, 
                        matched.dataset = match_statewide_PSM, 
                        full.dataset = dat.statewide, match.type = "PSM", 
                        model.adjust = "none")

r1

cov_to_adjust <- summary(m_statewide_PSM)$sum.matched[,3]
cov_to_adjust <- cov_to_adjust[names(cov_to_adjust) != "distance" & abs(cov_to_adjust) >= 0.1]

# adjust for covariates with SMD > 0.1
statewide_PSM_full_adjust <- glm(deaths ~ Any_bin + scale(bed.1000) + scale(percHisp) + scale(percWhite) +
                                   scale(days_since_fc) + scale(log(med_OOHH_value)) + scale(log(median_income)), 
                                 offset = log(tot_pop),
                                 family = "quasipoisson", data = match_statewide_PSM)

r1a <- extract_results.f(matched.model = statewide_PSM_full_adjust, 
                         matched.dataset = match_statewide_PSM, 
                         full.dataset = dat.statewide, match.type = "PSM", 
                         model.adjust = "adjusted further")

r1a

# PSM
m_ucmr_PSM <- matchit(
  Any5 ~ bed.1000 + percHisp + pop.dens.sqmi +
    percBlack + perc65plus + perchomeowner +
    log(med_OOHH_value) + percunderHSed + percWhite +
    log(median_income) + days_since_fc,
  data = dat.ucmr,
  method = "nearest",
  distance = "glm",
  exact = ~ census_region,
  caliper = 1,
  discard = "both",
  replace = FALSE,
  reestimate = TRUE)

match_ucmr_PSM <- match.data(m_ucmr_PSM, distance = "prop.score")
summary(m_ucmr_PSM)

# including variables with possible remaining imbalance also
ucmr_PSM <- glm(deaths ~ Any5, offset = log(tot_pop), 
                                  family = "quasipoisson", data = match_ucmr_PSM)

r2 <- extract_results.f(matched.model = ucmr_PSM, matched.dataset = match_ucmr_PSM, 
                        full.dataset = dat.ucmr, match.type = "PSM", 
                        model.adjust = "none")

r2

cov_to_adjust2 <- summary(m_ucmr_PSM)$sum.matched[,3]
cov_to_adjust2 <- cov_to_adjust2[names(cov_to_adjust2) != "distance" & abs(cov_to_adjust2) >= 0.1]

ucmr_PSM_full_adjust <- glm(deaths ~ Any5 + scale(bed.1000) + scale(percHisp) + 
                                          scale(pop.dens.sqmi) + scale(perchomeowner) + 
                                          scale(percWhite),
                                        offset = log(tot_pop),
                                        family = "quasipoisson", data = match_ucmr_PSM)

r2a <- extract_results.f(matched.model = ucmr_PSM_full_adjust, 
                         matched.dataset = match_ucmr_PSM, 
                        full.dataset = dat.ucmr, match.type = "PSM", 
                        model.adjust = "adjusted further")

r2a

all.objects <- ls()

all.results <- rbind(r1, r1a, r2, r2a)

all.results$dataset <- c(rep("Statewide sampling", 2), 
                         rep("UCMR 3", 2))

if (FALSE) {
  write.csv(all.results, "../results/matching_results.csv")
}

################################################################################
# love plots for visual assessment of covariate balance
v <- data.frame(old = c("days_since_fc", "log(med_OOHH_value)", 
                        "log(median_income)","bed.1000", "percHisp",
                        "census_region_New England", "pop.dens.sqmi",
                        "census_region_Middle Atlantic",
                        "census_region_Pacific", "percBlack", 
                        "census_region_South Atlantic", "percunderHSed",
                        "census_region_Mountain",
                        "census_region_East South Central", "census_region_East North Central",
                        "perc65plus", "perchomeowner", "percWhite"),
                new = c("Days since first case", "Median OOHH value (log)",
                        "Median income (log)", "Hospital beds",
                        "% Hispanic", "Region: New England", "Pop. density", 
                        "Region: Middle Atlantic", "Region: Pacific", "% Black",
                        "Region: South Atlantic", "% < HS education", 
                        "Region: Mountain", "Region: East South Central", 
                        "Region: East North Central", "% aged 65+",
                        "% homeowner", "% White")
)

love.plot(m_statewide_PSM, binary = "std",
          title = "Covariate balance for statewide data",
          thresholds = c(m = 0.1),
          shapes = c("square", "circle"),
          drop.distance = TRUE, position = "bottom",
          sample.names = c("Unmatched", "Matched"), 
          colors = MetBrewer::met.brewer("Egypt", 2),
          var.names = v,
          var.order = "unadjusted")

if (FALSE) {
  ggsave("../Figures/cov_balance_PSM_statewide.png", dpi = 300, width = 8, height = 8)
}

v <- data.frame(old = c("days_since_fc", "log(med_OOHH_value)", 
                        "log(median_income)","pop.dens.sqmi", "bed.1000",
                        "census_region_Middle Atlantic", "percHisp",
                        "census_region_New England",
                        "census_region_Pacific",
                        "census_region_South Atlantic",
                        "census_region_East North Central", "percBlack",
                        "census_region_Mountain", "perc65plus",
                        "percunderHSed", "census_region_East South Central",
                        "perchomeowner", "percWhite", 
                        "census_region_West North Central",
                        "census_region_West South Central"),
                new = c("Days since first case", "Median OOHH value (log)",
                        "Median income (log)", "Pop. density", "Hospital beds",
                        "Region: Middle Atlantic", "% Hispanic",
                        "Region: New England", "Region: Pacific", "Region: South Atlantic",
                        "Region: East North Central", "% Black", 
                        "Region: Mountain", "% aged 65+", "% < HS education",
                        "Region: East South Central", "% homeowner", "% White", 
                        "Region: West North Central", "Region: West South Central")
)

love.plot(m_ucmr_PSM, binary = "std",
          title = "Covariate balance for UCMR 3",
          thresholds = c(m = 0.1),
          drop.distance = TRUE, position = "bottom",
          shapes = c("square", "circle"),
          sample.names = c("Unmatched", "Matched"), 
          colors = MetBrewer::met.brewer("Egypt", 2),
          var.names = v,
          var.order = "unadjusted")

if (FALSE) {
  ggsave("../Figures/cov_balance_PSM_UCMR.png", dpi = 300, width = 8, height = 8)
}
