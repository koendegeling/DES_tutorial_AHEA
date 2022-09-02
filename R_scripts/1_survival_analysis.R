# This script belongs to the manuscript: Degeling K, Karnon J, Van de Ven M, Brennan A,
# Koffijberg H. Discrete event simulation in R using the 'simmer' package for health
# economic modeling: a tutorial and illustration in colon cancer. Applied Health 
# Economics and Health Policy.
#
# This is a supporting script that performs the parametric survival analysis of the 
# 'colon' data set from the 'survival' package, as well as United States life tables.
# to define the parameters for the time-to-event distributions that will be used in the
# discrete event simulation.
#
# In this and other scripts, we abbreviate to the different health states as follows:
#   - RF  recurrence free
#   - R   recurrence
#   - D   death
#
# Note that in the discrete event simulation, the competing events will be simulated
# using an event-specific distributions approach. This means that a time-to-event 
# distribution will be fitted for each competing event, which is censored for the
# competing event. More details on the data analysis and simulation for competing
# events in discrete event simulation can be found in: Degeling K et al. Comparing 
# strategies for modeling competing risks in discrete event simulations: a simulation
# study and illustration in colorectal cancer. Medical Decision Making. 2019. 
# - https://doi.org/10.1177/0272989X18814770"
# 
# The script has the following sections:
#   1. INITIALIZATION & DATA PREPATATION            loading packages and the data, and
#                                                   transforming the data in a useful 
#                                                   format to be analyzed
#   2. TIME-TO-RECURRENCE (TTR, RF -> R)            parametric models for TTR
#   3. BACKGROUND SURVIVAL (BS, RF -> D)            parametric models for BS
#   4. TIME FROM RECURRENCE TO DEATH (TTD, R -> D)  parametric models for TTD
#   5. EXPORT FITS                                  exporting the analysis outputs
#
# Version history:
# v1    28-10-2020  Koen Degeling     R 3.6.1   Initial version
# v2    20-07-2021  Koen Degeling     R 4.0.3   Updating comments




### 1. INITIALIZATION & DATA PREPARATION ----


# Uncomment to clear objects from Global Environment and clear the console
#rm(list = ls()); cat("\014");


# Packages
library(dplyr);         # version 1.0.6
library(tidyr);         # version 1.1.3
library(survival);      # version 3.2-7
library(flexsurv);      # version 2.0
library(flexsurvcure);  # version 1.2.0


# Data set
colon <- survival::colon;

# In the data set, each patient is represented by two rows, one for the event of 
# recurrence and one for the event of death. The variable "etype" indicates to which
# event the row correspond: 1 = recurrence and 2 = death. For each event, there is a
# "status" variable indicating whether the corresponding "time" corresponds to an 
# event observation or a censoring time. See the link below for more information:
# - https://cran.r-project.org/web/packages/survival/survival.pdf

# As mentioned before, the possible health states are defined to be:
# - RF  recurrence free
# - R   recurrence
# - D   death

# The following transitions are possible
# - RF -> R
# - RF -> D
# - R  -> D (only after recurrence)

# Hence, from the RF health state there are two competing risks: recurrence and death


# The first step is to transform the data into wide format, where the data from the
# etype, status and time variables is merged. Also, in the tutorial we only consider
# the two active treatment arms: rx = LEV and rx = Lev+5FU, so the observation strategy
# (rx = Obs) can be filtered out
colon_wide <- colon %>%  
  filter(rx != "Obs") %>% 
  pivot_wider(id_cols     = c(id, study, rx, sex, age, obstruct, perfor, adhere, nodes, differ, extent, surg, node4), 
              names_from  = etype, 
              values_from = c(time, status));


# For a traditional survival analysis in which each transition is modeled separately,
# two variables are defined for each transition: a censoring indicator "c" (1 = event, 
# 0 = censored) and the corresponding time "t".
df_colon <- colon_wide %>% 
  mutate(
    
    # Transition: RF -> R
    c_RF_R = status_1,
    t_RF_R = time_1,
    
    # Transition: RF -> D
    c_RF_D = case_when(
      status_1 == 1 ~ 0,
      status_2 == 0 ~ 0,
      status_1 == 0 & status_2 == 1 ~ 1,
      TRUE ~ NA_real_
    ),
    t_RF_D = if_else(status_1 == 1, time_1, time_2),
    
    # Transition: R -> D
    c_R_D  = if_else(status_1 == 0, NA_real_, status_2),
    t_R_D  = if_else(status_1 == 0, NA_real_, time_2 - time_1)
    
  );

# To tidy the data.frame up, the status and time variables are dropped
df_colon <- df_colon %>% 
  select(-contains("status"), -contains("time"));


# There are some "ties" (n=3) in the data with a recurrence being recorded at the same 
# day as death, which suggests an overall survival after progression of zero. Given
# that a zero is problematic for fitting time-to-even distributions these are manually
# corrected so that the event of recurrence happens 2 week (14 days) before the event 
# of death.
count(df_colon, (t_R_D == 0) & (c_R_D == 1));

df_colon <- df_colon %>% 
  mutate(
    correct_t_R_D = (t_R_D == 0) & (c_R_D == 1),
    t_RF_R = if_else(!is.na(correct_t_R_D) & correct_t_R_D, t_RF_R - 14, t_RF_R),
    t_R_D  = if_else(!is.na(correct_t_R_D) & correct_t_R_D, 14, t_R_D)
  );

head(df_colon %>% filter(correct_t_R_D));


# Changing time from days to years
days_in_year <- 365.25;
df_colon <- df_colon %>% 
  mutate(
    t_RF_R = t_RF_R / days_in_year,
    t_RF_D = t_RF_D / days_in_year,
    t_R_D  = t_R_D  / days_in_year
  );




### 2. TIME-TO-RECURRENCE (TTR, RF -> R) ----

# Quick inspection of S(t)
km_RF_R <- survfit(formula = Surv(t_RF_R, c_RF_R) ~ rx, data = df_colon);

plot(km_RF_R);

# As common for adjuvant treatment, there is a plateau at a certain point, as patients
# who have not progressed after 5 years are considered "cured" from a clinical point of
# view. Hence, a cure model is the most appropriate to use for this transition. These 
# can be fitted using the "flexsurvcure" package. For sake of ease, we parameterize the
# cure proportion theta and scale/rate/meanlog parameters only, not the shape/sdlog
# parameters.

# Similarly, for demonstration purposes, one covariate will be included in the survival
# models for RF -> R. To select the most influential covariate, a Cox proportional 
# hazards model is used: node4 is selected
coxph(formula = Surv(t_RF_R, c_RF_R) ~ rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + differ, 
      data    = df_colon);

# The probability of having more than 4 positive nodes in the study: 
# - Lev:      28.7%
# - Lev+5FU:  26.0%
# - Overall:  27.4% 
df_colon %>% 
  group_by(rx) %>% 
  summarize(p_node4 = mean(node4));

df_colon %>% 
  summarize(p_node4 = mean(node4));


# Fit distributions
fit_RF_R_exp      <- flexsurvcure(formula = Surv(t_RF_R, c_RF_R) ~ rx + node4 + rate(rx) + rate(node4), data = df_colon, dist = "exp");
fit_RF_R_gamma    <- flexsurvcure(formula = Surv(t_RF_R, c_RF_R) ~ rx + node4 + rate(rx) + rate(node4), data = df_colon, dist = "gamma");
fit_RF_R_gompertz <- flexsurvcure(formula = Surv(t_RF_R, c_RF_R) ~ rx + node4 + rate(rx) + rate(node4), data = df_colon, dist = "gompertz");
fit_RF_R_llogis   <- flexsurvcure(formula = Surv(t_RF_R, c_RF_R) ~ rx + node4 + scale(rx) + scale(node4), data = df_colon, dist = "llogis");
fit_RF_R_lnorm    <- flexsurvcure(formula = Surv(t_RF_R, c_RF_R) ~ rx + node4 + meanlog(rx) + meanlog(node4), data = df_colon, dist = "lnorm");
fit_RF_R_weibull  <- flexsurvcure(formula = Surv(t_RF_R, c_RF_R) ~ rx + node4 + scale(rx) + scale(node4), data = df_colon, dist = "weibull");

# Assess visual fit
# - Log-logistic looks best
par(mfrow = c(2, 3));
plot(fit_RF_R_exp, main = "exp");
plot(fit_RF_R_gamma, main = "gamma");
plot(fit_RF_R_gompertz, main = "gompertz");
plot(fit_RF_R_llogis, main = "llogis");
plot(fit_RF_R_lnorm, main = "lnorm");
plot(fit_RF_R_weibull, main = "weibull");

# Assess log-likelihood (higher/closer to zero is better)
# - Log-logistic best, but one parameter more than exponential
fit_RF_R_exp$loglik;
fit_RF_R_gamma$loglik;
fit_RF_R_gompertz$loglik;
fit_RF_R_llogis$loglik;
fit_RF_R_lnorm$loglik;
fit_RF_R_weibull$loglik;

# Assess AIC (lower is better)
# - Log-logistic is best
AIC(fit_RF_R_exp);
AIC(fit_RF_R_gamma);
AIC(fit_RF_R_gompertz);
AIC(fit_RF_R_llogis);
AIC(fit_RF_R_lnorm);
AIC(fit_RF_R_weibull);

# Assess BIC (lower is better)
# - Log-logistic is best
BIC(fit_RF_R_exp);
BIC(fit_RF_R_gamma);
BIC(fit_RF_R_gompertz);
BIC(fit_RF_R_llogis);
BIC(fit_RF_R_lnorm);
BIC(fit_RF_R_weibull);

# Log-logistic selected based on best likelihood-based measures.
fit_RF_R <- fit_RF_R_llogis;




### 3. BACKGROUND SURVIVAL (BS, RF -> D) ----

# Quick inspection of S(t)
km_RF_D <- survfit(formula = Surv(t_RF_D, c_RF_D) ~ rx, data = df_colon);

par(mfrow = c(1, 1)); plot(km_RF_D);

# Univariable
km_RF_D_uni <- survfit(formula = Surv(t_RF_D, c_RF_D) ~ 1, data = df_colon);

plot(km_RF_D_uni);

# No difference between arm (as expected for background mortality) and limited follow-
# up results in low number of events. Assuming a Gompertz distribution makes more sense
# than fitting separate distributions and picking based on likelihood-based measures.

fit_RF_D <- flexsurvreg(Surv(t_RF_D, c_RF_D) ~ 1, data = df_colon, dist = "Gompertz");

q_sim <- seq(0, 60, 0.1);
plot(km_RF_D_uni, xlim = c(0, max(q_sim)));
lines(x = q_sim, 
      y = pgompertz(q = q_sim, shape = fit_RF_D$coefficients["shape"], rate = exp(fit_RF_D$coefficients["rate"]), lower.tail = F),
      col = "red");

# Survival is still too high given that the mean age at randomization in the trial was
# 60 years. Hence, another approach is needed: using life tables.

# Lifetables for the US were obtained from the National Vital Statistics Reports - 
# United States Life Tables, 2017, by Elizabeth Arias, Ph.D., and Jiaquan Xu, M.D., 
# Division of Vital Statistics - Table B on page 4
# - webpage: https://www.cdc.gov/nchs/products/life_tables.htm#life
# - report:  https://www.cdc.gov/nchs/data/nvsr/nvsr68/nvsr68_07-508.pdf
df_US_lifetable_2017 <- data.frame(
  age = c(60,    65,    70,    75,    80,    85,    90,    95,   100),
  n   = c(88226, 83696, 77697, 69418, 57839, 42382, 24560, 9361, 1894),
  t   = c(0,     5,     10,    15,    20,    25,    30,    35,   40)
) %>% mutate(
  p = n / n[1]
);

# A Gompertz is fitted by non-linear least squares regression, using the parameters of
# the previously fitted distribution as start values. The main downside of this 
# approach is that parameter uncertainty cannot be reflected as the parameter estimates
# are not estimated by maximum likelihood. However, the uncertainty in the life tables 
# is very low, because of the large number of individuals on which they are based, so 
# this is not relevant in this case.
fit_RF_D_lifetable <- nls(formula = p ~ pgompertz(q = t, shape = x_shape, rate = exp(x_rate), lower.tail = FALSE),
                          start   = list(x_shape = fit_RF_D$coefficients["shape"],
                                         x_rate  = fit_RF_D$coefficients["rate"]),
                          data    = df_US_lifetable_2017);

# As the plot shows, this results in a much more realistic curve
q_sim <- seq(0, 60, 0.1);
plot(km_RF_D_uni, xlim = c(0, max(q_sim)));
lines(x = df_US_lifetable_2017$t, y = df_US_lifetable_2017$p, col = "blue")
lines(x = q_sim, 
      y = pgompertz(q = q_sim, shape = coef(fit_RF_D_lifetable)["x_shape"], rate = exp(coef(fit_RF_D_lifetable)["x_rate"]), lower.tail = F),
      col = "red",
      lty = 2);




### 4. TIME FROM RECURRENCE TO DEATH (TTD, R -> D) ----

# Note that in the discrete event simulation, we will use TTD as time from recurrence
# to death due to cancer, which is a competing event to death due to other causes, for
# which we estimated the Gompertz distribution in the previous section. Technically, 
# this is not true, because we cannot single out cancer-related death from potential
# non-cancer-related deaths in the dataset. However, given that this tutorial is about
# how to implement discrete event simulation in R, rather than on survival analysis,
# this solution is appropriate. Also, patients with a colorectal cancer recurrence are
# most likely to die from colorectal cancer-related causes, so it is a plausible 
# assumption to interpret the data this way.

# Quick inspection of S(t)
km_R_D <- survfit(formula = Surv(t_R_D, c_R_D) ~ rx, data = df_colon);

par(mfrow = c(1, 1)); plot(km_R_D, col = c("black", "blue", "green"));
km_R_D$strata;


# Again, common trend. Use of standard parametric models makes most sense. Might be 
# unexpected that the most effective regimen with regard to TTR results in lowest post-
# recurrence OS (Lev+5FU), put those probably represent the most aggressive cancers. 
# but proportional models seems problematic, so parameterizing all parameters of the
# time-to-event distributions.
# - note that there is one individual without an event and a time of zero, which causes
#   problems and, hence, this observation is filtered out
df_colon %>% 
  count(t_R_D == 0)

fit_R_D_exp      <- flexsurvreg(formula = Surv(t_R_D, c_R_D) ~ rx, data = filter(df_colon, t_R_D > 0), dist = "exp");
fit_R_D_gamma    <- flexsurvreg(formula = Surv(t_R_D, c_R_D) ~ rx + shape(rx), data = filter(df_colon, t_R_D > 0), dist = "gamma");
fit_R_D_gompertz <- flexsurvreg(formula = Surv(t_R_D, c_R_D) ~ rx + shape(rx), data = filter(df_colon, t_R_D > 0), dist = "gompertz");
fit_R_D_llogis   <- flexsurvreg(formula = Surv(t_R_D, c_R_D) ~ rx + shape(rx), data = filter(df_colon, t_R_D > 0), dist = "llogis");
fit_R_D_lnorm    <- flexsurvreg(formula = Surv(t_R_D, c_R_D) ~ rx + sdlog(rx), data = filter(df_colon, t_R_D > 0), dist = "lnorm");
fit_R_D_weibull  <- flexsurvreg(formula = Surv(t_R_D, c_R_D) ~ rx + shape(rx), data = filter(df_colon, t_R_D > 0), dist = "weibull");

# Fits are not great, though Log-logistic seems best
par(mfrow = c(2, 3));
plot(fit_R_D_exp, main = "exp");
plot(fit_R_D_gamma, main = "gamma");
plot(fit_R_D_gompertz, main = "gompertz");
plot(fit_R_D_llogis, main = "llogis");
plot(fit_R_D_lnorm, main = "lnorm");
plot(fit_R_D_weibull, main = "weibull");

# Log-logistic best based on log-likelihood
fit_R_D_exp$loglik;
fit_R_D_gamma$loglik;
fit_R_D_gompertz$loglik;
fit_R_D_llogis$loglik;
fit_R_D_lnorm$loglik;
fit_R_D_weibull$loglik;

# Log-logistic best based on AIC
AIC(fit_R_D_exp);
AIC(fit_R_D_gamma);
AIC(fit_R_D_gompertz);
AIC(fit_R_D_llogis);
AIC(fit_R_D_lnorm);
AIC(fit_R_D_weibull);

# Log-logistic best based on BIC
BIC(fit_R_D_exp);
BIC(fit_R_D_gamma);
BIC(fit_R_D_gompertz);
BIC(fit_R_D_llogis);
BIC(fit_R_D_lnorm);
BIC(fit_R_D_weibull);        

# Log-logistic selected based on visual fit and likelihood-based measures.
fit_R_D <- fit_R_D_llogis;




### 5. EXPORT FITS ----

# The parameters are exported to a CSV file with strategies indicated by numbers:
#   1 Lev
#   2 Lev+5FU
# And the number of nodes indicated by the corresponding risk of recurrence:
#   low   4 or less positive nodes
#   high  more than 4 positive nodes

# Functions to transform coefficients to parameters on real scale:
# - getProb     to transform the theta parameter used in the "flexsurvcure" package,
#               which is on logit scale, to a cure probability/proportion
# - getExpPar   to transform parameters on log scale to real scale
getProb   <- function(x) unname(exp(sum(x)) / (1 + exp(sum(x))));
getExpPar <- function(x) unname(exp(sum(x)));

# First combine in one big list
ls_survival_models <- list(
  
  # RF -> R (for the llogis both the shape and scale need to be exponentiated)
  d_RF_R_1_low = list(dist  = "llogis cure", 
                      cure  = getProb(fit_RF_R$coefficients["theta"]),
                      shape = getExpPar(fit_RF_R$coefficients["shape"]),
                      scale = getExpPar(fit_RF_R$coefficients["scale"])),
  d_RF_R_1_high = list(dist  = "llogis cure", 
                       cure  = getProb(fit_RF_R$coefficients[c("theta", "node4")]),
                       shape = getExpPar(fit_RF_R$coefficients[c("shape")]),
                       scale = getExpPar(fit_RF_R$coefficients[c("scale", "scale(node4)")])),
  d_RF_R_2_low = list(dist  = "llogis cure", 
                      cure  = getProb(fit_RF_R$coefficients[c("theta", "rxLev+5FU")]),
                      shape = getExpPar(fit_RF_R$coefficients[c("shape")]),
                      scale = getExpPar(fit_RF_R$coefficients[c("scale", "scale(rxLev+5FU)")])),
  d_RF_R_2_high = list(dist  = "llogis cure", 
                       cure  = getProb(fit_RF_R$coefficients[c("theta", "rxLev+5FU", "node4")]),
                       shape = getExpPar(fit_RF_R$coefficients[c("shape")]),
                       scale = getExpPar(fit_RF_R$coefficients[c("scale", "scale(rxLev+5FU)", "scale(node4)")])),
  
  # RF -> D (for the gompertz only the rate needs to be exponentiated)
  d_RF_D = list(dist  = "gompertz",
                shape = unname(coef(fit_RF_D_lifetable)["x_shape"]),
                rate  = getExpPar(coef(fit_RF_D_lifetable)["x_rate"])),
  
  # R -> D (for the llogis both the shape and scale need to be exponentiated)
  d_R_D_1 = list(dist  = "llogis",
                 shape = getExpPar(fit_R_D$coefficients["shape"]),
                 scale = getExpPar(fit_R_D$coefficients["scale"])),
  d_R_D_2 = list(dist  = "llogis",
                 shape = getExpPar(fit_R_D$coefficients[c("shape", "shape(rxLev+5FU)")]),
                 scale = getExpPar(fit_R_D$coefficients[c("scale", "rxLev+5FU")]))
  
);

# Transform to a data.frame
df_survival_models <- as.data.frame(ls_survival_models);

# Correct column names by replacing "." by "_"
colnames(df_survival_models) <- gsub(pattern = "\\.", replacement = "_", x = colnames(df_survival_models));

# Export data.frame
write.csv(x = df_survival_models, file = "data/survival_models.csv", row.names = FALSE);

# Export the survival analysis objects to sample parameter values for the probabilistic 
# analysis
save(fit_RF_R, fit_R_D, file = "data/survival_models.RDS");



