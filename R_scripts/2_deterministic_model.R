# This script belongs to the manuscript: Degeling K, Karnon J, Van de Ven M, Brennan A,
# Koffijberg H. Discrete event simulation in R using the 'simmer' package for health
# economic modeling: a tutorial and illustration in colon cancer. Applied Health 
# Economics and Health Policy.
# 
# This is the main script of the tutorial, demonstrating the implementation and 
# deterministic analysis of the discrete event simulation (DES). Note that the aim of 
# this script is to provide a gentle introduction to implementing DES using the 
# "simmer" package and, hence, the code is written to be as straightforward to 
# interpret as possible, which may not necessarily represent the most efficient coding
# approach.
#
# Please refer to the manuscript for further information about the case study and the 
# model parameters, as well as to its Appendix for an extensive introduction to the
# function from the "simmer" package that are used.
#
# One important aspect to highlight for those who are not interested in understanding
# the details of the "simmer" package and its function, but who do want to obtain a
# high-level understanding of this script, is that DES models that are implemented 
# using this package exist of two main object: the trajectory and the simulation
# environment. Trajectories define the model structure, i.e. the process to which the
# individuals flow. Simulation environments define the number of individuals to be
# simulated through a certain trajectory and tracks the state of the simulation and the
# individuals over the simulation time.
#
# When using the "simmer" package, the functions, parameters and potential other 
# objects that are referred to in the trajectory, have to be defined before defining 
# the trajectory. Therefore, this script is structured as follows:
#   1. INITIALIZATION           loading packages and survival analysis results
#   2. MODEL PARAMETERS         defining the parameters used in the model
#   3. SUPPORTING FUNCTIONS     defining custom functions for a  tidy implementation
#   4. MODEL STRUCTURE          definition of the trajectory
#   5. RUN SIMULATION           definition and running the simulation
#   6. STABILITY OF OUTCOMES    assessing the number of individuals to be simulated
#
# Some further information that may be helpful:
# - Time is defined in years
# - The snake_case is used to identify parameters and objects (e.g., d_TTD_shape_1 and
#   fn_adjuvant_cycle) and the CamelCase is used to identify attributes/variables in
#   the simulation (e.g. TreatmentArm and AdjuvantCycles). This allows to easily 
#   distinguish between what are model parameter or functions, and what are attributes.
#
# Version history:
# v1    28-10-2020  Koen Degeling     R 3.6.1   Initial version
# v2    20-07-2021  Koen Degeling     R 4.0.3   Updating comments




### 1. INITIALIZATION ----

# Uncomment to clear objects from Global Environment and clear the console
#rm(list = ls()); cat("\014");

# Required packages
library(simmer);        # version 4.4.2
library(simmer.plot);   # version 0.1.16
library(flexsurv);      # version 2.0 - for sampling from distributions
library(data.table);    # version 1.14.0 - for summarizing the simulation output

# Loading a data.frame that contains the parameters estimated in the survival analysis
df_survival_models <- read.csv(file = "data/survival_models.csv");




### 2. MODEL PARAMETERS ----

# The first part of the parameter names identifies the type of parameter:
#   c_    cost
#   d_    distribution parameter
#   dr_   discount rate
#   n_    number/count
#   p_    probabilty
#   t_    time/duration
#   u_    health state utility

# Where applicable, the last part of the parameter names identifies the strategy to 
# which the parameter applies:
#   _1    control strategy: Lev
#   _2    experimental strategy: Lev + 5FU

# Where applicable, the last part of the parameter names identifies the risk group
# to which the parameter is applicable:
#   low    low risk
#   high   high risk


## GENERAL PARAMETERS / CHARACTERISTICS / BACKGROUND SURVIVAL

# Discount rates
dr_costs  <- 0.04;
dr_health <- 0.015;

# Probability of being high risk, defined by whether more than 4 lymph nodes were found
# to be positive (high risk) or not (low risk), which is used as a covariate for 
# time-to-recurrence (TTR). The proportion of high-risk patients was 27.4% in the 
# original study, but here were are interested in a different population.
p_highrisk <- 0.6;


# Gompertz distribution for death due to other causes
d_BS_shape <- df_survival_models$d_RF_D_shape;
d_BS_rate  <- df_survival_models$d_RF_D_rate;


## TIME-TO-RECURRENCE / ADJUVANT TREATMENT

# Log-logistic cure distribution for time-to-recurrence (TTR)
# - conditional on treatment arm and risk group
d_TTR_cure_1_low  <- df_survival_models$d_RF_R_1_low_cure;
d_TTR_cure_1_high <- df_survival_models$d_RF_R_1_high_cure;
d_TTR_cure_2_low  <- df_survival_models$d_RF_R_2_low_cure;
d_TTR_cure_2_high <- df_survival_models$d_RF_R_2_high_cure;

d_TTR_shape_1_low  <- df_survival_models$d_RF_R_1_low_shape;
d_TTR_shape_1_high <- df_survival_models$d_RF_R_1_high_shape;
d_TTR_shape_2_low  <- df_survival_models$d_RF_R_2_low_shape;
d_TTR_shape_2_high <- df_survival_models$d_RF_R_2_high_shape;

d_TTR_scale_1_low  <- df_survival_models$d_RF_R_1_low_scale;
d_TTR_scale_1_high <- df_survival_models$d_RF_R_1_high_scale;
d_TTR_scale_2_low  <- df_survival_models$d_RF_R_2_low_scale;
d_TTR_scale_2_high <- df_survival_models$d_RF_R_2_high_scale;

# Adjuvant treatment costs
c_adjuvant_cycle_1 <- 5000;
c_adjuvant_cycle_2 <- 15000;

# Utility during adjuvant treatment
u_adjuvant <- 0.70;

# Other adjuvant treatment parameters
t_adjuvant_cycle      <- 3/52; # in years
n_max_adjuvant_cycles <- 10;

# Probability of toxicities during a cycle
# - conditional on treatment arm
p_tox_1 <- 0.20;
p_tox_2 <- 0.40;

# Costs of toxicities
c_tox <- 2000;

# Disutility of toxicities
u_dis_tox <- 0.10;


## DISEASE MONITORING

# Cost of a monitoring cycle
c_monitor <- 1000;

# Utility while free of recurrence
u_diseasefree <- 0.80;

# Other disease monitoring parameters
t_monitor_cycle      <- 1; # in years
n_max_monitor_cycles <- 5;


## RECURRENCE OF DISEASE

# Log-logistic distribution for the time-to-death due to cancer after recurrence (TTD)
# - conditional on treatment arm
d_TTD_shape_1 <- df_survival_models$d_R_D_1_shape;
d_TTD_shape_2 <- df_survival_models$d_R_D_2_shape;

d_TTD_scale_1 <- df_survival_models$d_R_D_1_scale;
d_TTD_scale_2 <- df_survival_models$d_R_D_2_scale;

# Cost of treatment for advanced disease
c_advanced <- 40000;

# Utility during advanced disease
u_advanced <- 0.60;




### 3. SUPPORTING FUNCTIONS ----

# A range of functions are defined to support the implementation and analysis of the
# model. Generally speaking, all custom function are identified their name starting 
# with "fn_", with a dew exceptions.

# Three types of supporting functions are defined:
#   1) General functions:     ... not specific to this model
#   2) Short-cut functions:   ... not specific to this model, used to tidy up the code
#   3) Modeling functions:    ... too complex/large to nicely include in the structure


## GENERAL FUNCTIONS 

# Functions for discounting
fn_discount_costs <- function(amount, t) amount / (1 + dr_costs) ^ t;
fn_discount_QALYs <- function(utility, t_start, t_duration) utility / (1 + dr_health) ^ t_start / log(1 + dr_health) * (1 - 1 / (1 + dr_health) ^ t_duration);

# Function that summarizes the recorded attribute values obtained from the simmer 
# environment using the get_mon_attributes() function. This is necessary because the 
# get_mon_attributes() function returns a data.frame with multiple rows per individual
# and potentially multiple rows per individual per attribute if the value of that 
# attribute changes over time. For the outcomes of the analysis, we are interested in
# the final value of the attributes, which is what this function conveniently obtains.
# - if keys (i.e. the attributes to summarize) are not specified, all patient-level
#   attributes will be summarized
# - it uses the data.table package because the data.frame returned by the
#   get_mon_attributes() function can be pretty large
fn_summarize <- function(df_sim, keys = NULL) {
  
  df <- as.data.table(df_sim)[name != ""]; # transform to data.table and remove any global attributes
  
  if(is.null(keys)) keys <- unique(df$key); # if no keys specified, use all keys
  
  df <- df[key %in% keys]; # filter to include only the keys of interest
  setorder(df, name, time); # order by individual and simulation time
  df <- df[, .(value = value[.N]), by = list(name, key)]; # group by individual and key and select the last value
  df <- dcast(df, name~key, value.var = "value"); # transform to wide format
  setcolorder(df, c("name", keys)); # order the columns
  
  return(df);
  
}

# Function to sample from a log-logistic cure distribution
# - it is vectorized, though this is slower than "if() ... else ..." when n = 1, so if
#   you want to speed up the simmer simulation (which is not vectorized), you can do so
#   by replacing "ifelse(..., ..., ...)" by "if(...) {...} else {...}"
# - you could also use one rather than two random numbers and sample using the qllogis
#   function for those without cure after rescaling the random number
rllogiscure <- function(n = 1, cure, shape, scale) ifelse(runif(n) < cure, Inf, rllogis(n = n, shape = shape, scale = scale));


## SHORT-CUT FUNCTIONS

# To tidy up the code, some short-cut functions are defined to prevent a lot of
# repetition in the model structure. The main thing that these short cuts work around
# is the defining the simulation environment through the .env argument. The simulation
# environment is defined as "sim" and included as default value for the corresponding
# argument in the short-cut functions. Furthermore, the names are shortened a bit .
set_attr <- set_attribute;
get_attr <- function(keys, .env = sim) get_attribute(.env = .env, keys = keys);
now_     <- function(.env = sim) now(.env = .env);


## MODELLING FUNCTIONS

fn_initialisation <- function() {
  
  # This function is used to provide the initial values for a range of attributes. It
  # has no specific input arguments, but uses a set of global parameters/variables.
  
  # Output:
  # [1] TreatmentArm
  # [2] HighRisk
  # [3] AdjuvantCycles
  # [4] MonitorCycles
  # [5] Toxicities
  # [6] dCosts
  # [7] dQALYs
  # [8] BS (background survival)
  # [9] TTR (time-to-recurrence)

  # The treatment arm is set based on the parameter treatment_arm that is specified 
  # before running the simulation.
  TreatmentArm <- treatment_arm;
  
  # Sampling whether the individual is at high risk of recurrence (1) or not (0).
  HighRisk <- as.integer(runif(1) < p_highrisk);

  # Initializing counters
  AdjuvantCycles <- 0;
  MonitorCycles <- 0;
  Toxicities <- 0;
  dCosts <- 0;
  dQALYs <- 0;
  
  # Sampling BS
  BS <- rgompertz(n = 1, shape = d_BS_shape, rate = d_BS_rate);
  
  # Sampling TTR conditional on TreatmentArm and HighRisk
  if(TreatmentArm == 1 & HighRisk == 0) {
    TTR <- rllogiscure(cure   = d_TTR_cure_1_low, 
                       shape  = d_TTR_shape_1_low, 
                       scale  = d_TTR_scale_1_low);
  } else if(TreatmentArm == 1 & HighRisk == 1) {
    TTR <- rllogiscure(cure   = d_TTR_cure_1_high, 
                       shape  = d_TTR_shape_1_high, 
                       scale  = d_TTR_scale_1_high);
  } else if(TreatmentArm == 2 & HighRisk == 0) {
    TTR <- rllogiscure(cure   = d_TTR_cure_2_low, 
                       shape  = d_TTR_shape_2_low, 
                       scale  = d_TTR_scale_2_low);
  } else if(TreatmentArm == 2 & HighRisk == 1) {
    TTR <- rllogiscure(cure   = d_TTR_cure_2_high, 
                       shape  = d_TTR_shape_2_high, 
                       scale  = d_TTR_scale_2_high);
  } else { 
    stop("This should not happen in selecting the TTR distribution!");
  }
  
  # Define the output object in the required order
  output <- c(TreatmentArm, HighRisk, AdjuvantCycles, MonitorCycles, Toxicities, dCosts, dQALYs, BS, TTR);
  
  return(output)
  
}

fn_adjuvant_cycle <- function(attrs, t_now) {
  
  # This function is used to determine what will happen to the individual in the 
  # upcoming adjuvant treatment cycle. It uses a range of individual-level attributes  
  # provided through the attrs argument, as well as the current simulation time through
  # provided the t_now argument.
  
  # Example use: 
  # fn_adjuvant_cycle(attrs = get_attr(keys = c("BS", "TTR")),  t_now = now_())
  
  # Input:
  # - attrs
  #   [1] BS
  #   [2] TTR
  # - t_now
  BS  <- attrs[1];
  TTR <- attrs[2];

  # Output:
  # [1] AdjuvantCycleEvent:
  #       1 = death during cycle
  #       2 = recurrence during cycle
  #       3 = no recurrence or death during cycle 
  # [2] AdjuvantCycleTime
  
  # Determining when a full cycle would end
  t_full_cycle_end <- t_now + t_adjuvant_cycle;
  
  # Determine the event and corresponding time-to-event
  # 1) If the patient dies during the cycle
  if((BS <= t_full_cycle_end) & (BS <= TTR)) {
    AdjuvantCycleEvent <- 1;
    AdjuvantCycleTime  <- BS - t_now;
  
  # 2) If the cancer recurs during the cycle
  } else if((TTR <= t_full_cycle_end) & (TTR < BS)) {
    AdjuvantCycleEvent <- 2;
    AdjuvantCycleTime  <- TTR - t_now;
  
  # 3) If the patient completes the cycle
  } else {
    AdjuvantCycleEvent <- 3;
    AdjuvantCycleTime  <- t_adjuvant_cycle;
  }
  
  # Define the output object in the required order
  output <- c(AdjuvantCycleEvent, AdjuvantCycleTime);
  
  return(output);
  
}

fn_adjuvant_impact <- function(attrs, t_now) {
  
  # This function is used to update the value of several counter attributes. It uses a
  # range of individual-level attributes provided through the attrs argument, as well 
  # as the current simulation time provided through the t_now argument.
  
  # NOTE! The outputs of this function are increments, i.e. they are added to the
  # respective counter attributes rather than replacing them. This is specified
  # through the 'mod = "+"' argument to the set_attr() function.
  
  # Example use:
  # fn_adjuvant_impact(attrs = get_attr(keys = c("TreatmentArm", "AdjuvantCycleTime")),
  #                    t_now = now_())
  
  # Input:
  # - attrs    (vector of attribute(s) previously set within the trajectory)
  #   [1] TreatmentArm
  #   [2] AdjuvantCycleTime
  # - t_now
  TreatmentArm      <- attrs[1];
  AdjuvantCycleTime <- attrs[2];
  
  # Output:
  # [1] AdjuvantCycles
  # [2] Toxicities
  # [3] dCosts
  # [4] dQALYs
  
  # The increment for the number of cycles is 1
  AdjuvantCycles <- 1;
  
  # Determine whether a toxicity will occur based on TreatmentArm
  p_toxicity <- if(TreatmentArm == 1) p_tox_1 else p_tox_2;
  Toxicities <- if(runif(1) < p_toxicity) 1 else 0;
  
  # Determine the increment in costs based on TreatmentArm and Toxicities
  c_adjuvant_cycle <- if(TreatmentArm == 1) c_adjuvant_cycle_1 else c_adjuvant_cycle_2;
  dCosts  <- fn_discount_costs(amount = c_adjuvant_cycle + Toxicities*c_tox,
                               t      = t_now);
  
  # Determine the increment in QALYs based on TreatmentArm and Toxicities
  dQALYs <- fn_discount_QALYs(utility    = u_adjuvant - Toxicities*u_dis_tox,
                              t_start    = t_now,
                              t_duration = AdjuvantCycleTime);
  
  # Define the output object in the required order
  output <- c(AdjuvantCycles, Toxicities, dCosts, dQALYs);
  
  return(output);
  
}

fn_monitor_cycle <- function(attrs, t_now) {
  
  # This function is used to determine what will happen to the individual in the 
  # upcoming disease monitoring cycle. It uses a range of individual-level attribute  
  # provided through the attrs argument, as well as the current simulation time through
  # provided the t_now argument.
  
  # Example use:
  # fn_monitor_cycle(attrs = get_attr(keys = c("MonitorCycles", "BS", "TTR")),
  #                  t_now = now_())
  
  # Input:
  # - attrs    (vector of attribute(s) previously set within the trajectory)
  #   [1] MonitorCycles
  #   [2] BS
  #   [2] TTR
  # - t_now
  MonitorCycles <- attrs[1];
  BS  <- attrs[2];
  TTR <- attrs[3];
  
  # Output:
  # [1] MonitorCycleEvent:
  #       1 = death during cycle
  #       2 = recurrence during cycle
  #       3 = no recurrence or death during cycle 
  # [2] MonitorCycleTime
  
  # Determining when a full cycle would end
  t_full_cycle_end <- t_now + t_monitor_cycle;
  
  # Determine the event and corresponding time-to-event
  # 1) If the patient dies during the cycle
  if((BS <= t_full_cycle_end) & (BS <= TTR)) {
    MonitorCycleEvent <- 1;
    MonitorCycleTime  <- BS - t_now;
  
  # 2) If the cancer recurs during the cycle
  } else if((TTR <= t_full_cycle_end) & (TTR < BS)) {
    MonitorCycleEvent <- 2;
    MonitorCycleTime  <- TTR - t_now;
  
  # 3) If the patient completes the cycle
  } else {
    MonitorCycleEvent <- 3;
    MonitorCycleTime  <- t_monitor_cycle;
  }
  
  # Define the output object in the required order
  output <- c(MonitorCycleEvent, MonitorCycleTime);
  
  return(output);
  
}

fn_monitor_impact <- function(attrs, t_now) {
  
  # This function is used to update the value of several counter attributes. It uses a
  # range of individual-level attributes provided through the attrs argument, as well 
  # as the current simulation time provided through the t_now argument.
  
  # NOTE! The outputs of this function are increments, i.e. they are added to the
  # respective counter attributes rather than replacing them. This is specified
  # through the 'mod = "+"' argument to the set_attr() function.
  
  # Example use:
  # fn_monitor_impact(attrs = get_attr(keys = c("MonitorCycleTime")), t_now = now_())
  
  # Input:
  # - attrs    (vector of attribute(s) previously set within the trajectory)
  #   [1] MonitorCycleTime
  # - t_now
  MonitorCycleTime <- attrs[1];
  
  # Output:
  # [1] MonitorCycles
  # [2] dCosts
  # [3] dQALYs
  
  # The increment for the number of cycles is 1
  MonitorCycles <- 1;

  # Determine the increment in costs
  dCosts  <- fn_discount_costs(amount = c_monitor,
                               t      = t_now + MonitorCycleTime);
  
  # Determine the increment in QALYs
  dQALYs <- fn_discount_QALYs(utility    = u_diseasefree,
                              t_start    = t_now,
                              t_duration = MonitorCycleTime);
  
  # Define the output object in the required order
  output <- c(MonitorCycles, dCosts, dQALYs);
  
  return(output);
  
}

fn_advanced_time <- function(attrs) {
  
  # This function determine the time the individual will spend on advanced treatment
  # after a recurrence was experienced. It uses a range of individual-level attributes 
  # provided through the attrs argument.
  
  # Example use:
  # fn_advanced_time(attrs = get_attr(keys = c("TreatmentArm", "BS")))
  
  # Input:
  # - attrs    (vector of attribute(s) previously set within the trajectory)
  #   [1] TreatmentArm
  #   [2] BS
  TreatmentArm <- attrs[1];
  BS <- attrs[2];
  
  # Output:
  # [1] TTD

  # Sample TTD based on TreatmentArm
  if(TreatmentArm == 1) {
    TTD <- rllogis(n     = 1,
                   shape = d_TTD_shape_1,
                   scale = d_TTD_scale_1);
  } else if(TreatmentArm == 2) {
    TTD <- rllogis(n     = 1,
                   shape = d_TTD_shape_2,
                   scale = d_TTD_scale_2);
  } else { 
    stop("This should not happen in selecting the TTD distribution!");
  }
  
  # Check whether TTD is lower than BS, if not correct to BS.
  TTD <- min(BS, TTD);
  
  # Define the output object in the required order
  output <- c(TTD);
  
  return(output);
  
}

fn_advanced_impact <- function(attrs, t_now) {
  
  # This function is used to update the value of several counter attributes. It uses a
  # range of individual-level attributes provided through the attrs argument, as well
  # as the current simulation time provided through the t_now argument.
  
  # NOTE! The outputs of this function are increments, i.e. they are added to the
  # respective counter attributes rather than replacing them. This is specified
  # through the 'mod = "+"' argument to the set_attr() function.
  
  # Example use:
  # fn_advanced_impact(attrs = get_attr(keys  = "TTD"), t_now = now_())
  
  # Input:
  # - attrs    (vector of attribute(s) previously set within the trajectory)
  #   [1] TTD
  # - t_now
  TTD <- attrs[1];
  
  # Output:
  # [1] dCosts
  # [2] dQALYs
  
  # Determine the increment in costs
  dCosts  <- fn_discount_costs(amount = c_advanced,
                               t      = t_now);
  
  # Determine the increment in QALYs
  dQALYs <- fn_discount_QALYs(utility    = u_advanced,
                              t_start    = t_now,
                              t_duration = TTD);
  
  # Define the output object in the required order
  output <- c(dCosts, dQALYs);
  
  return(output);
  
}




### 4. MODEL STRUCTURE ----

# The main model structure is defined by the trajectory object traj_main. Given that
# there are some parts of the code that are repetitive, including recording survival
# when the individual dies and treatment of advanced disease, these sections of code 
# are defined in sub-trajectories: traj_death and traj_recurrence, respectively. These
# need to be available from the global R environment when defining traj_main and, hence,
# they are defined first even though that may seem counterintuitive from a disease
# pathway/chronological point of view.
#
# Where applicable, the numbering in the comments corresponds to the numbering in the
# pseudocode flowchart presented in the manuscript and the manuscript sub-sections.


# Sub-trajectory to record survival upon death
traj_death <- trajectory(name = "traj_death") %>% 
  
  # 17) Record survival
  set_attr(keys = "OS", values = function() now_());


# Sub-trajectory for treatment of recurrence/advanced disease
traj_recurrence <- trajectory(name = "traj_recurrence") %>% 
  
  ## TREATMENT OF ADVANCED DISEASE ##
  
  # 14) Record the time until death
  set_attr(keys   = "TTD",
           values = function() fn_advanced_time(attrs = get_attr(keys = c("TreatmentArm", "BS")))) %>% 
  
  # 15) Update the costs and QALYs
  set_attr(keys   = c("dCosts", "dQALYs"),
           values = function() fn_advanced_impact(attrs = get_attr(keys  = "TTD"),
                                                  t_now = now_()),
           mod    = "+") %>% 
  
  # 16) Delay until death
  timeout_from_attribute(key = "TTD") %>% 
  
  # Death: go to traj_death sub-trajectory
  join(traj_death);


# Main trajectory
traj_main <- trajectory(name = "traj_main") %>%
  
  ## INITIALIZATION ##
  
  # 1) Initialize the individual-level attributes
  set_attr(keys   = c("TreatmentArm", "HighRisk", "AdjuvantCycles", "MonitorCycles", "Toxicities", "dCosts", "dQALYs", "BS", "TTR"),
           values = function() fn_initialisation()) %>% 
  
  
  ## ADJUVANT TREATMENT ##
  
  # 2) Record the event and duration for the treatment cycle
  set_attr(keys   = c("AdjuvantCycleEvent", "AdjuvantCycleTime"),
           values = function() fn_adjuvant_cycle(attrs = get_attr(keys = c("BS", "TTR")),
                                                 t_now = now_())) %>% 
  
  # 3) Update the number of cycles, toxicities, costs, and QALYs
  set_attr(keys   = c("AdjuvantCycles", "Toxicities", "dCosts", "dQALYs"),
           values = function() fn_adjuvant_impact(attrs = get_attr(keys = c("TreatmentArm", "AdjuvantCycleTime")),
                                                  t_now = now_()),
           mod    = "+") %>% 
  
  # 4) Delay for the duration of the treatment cycle
  timeout_from_attribute(key = "AdjuvantCycleTime") %>% 
  
  # 5) Check what will happen next based on the previously set AdjuvantCycleEvent:
  #     1 = death during cycle
  #     2 = recurrence during cycle
  #     3 = no recurrence or death during cycle 
  ## 5.1) Did the individual die?
  branch(option = function() (get_attr(keys = "AdjuvantCycleEvent") == 1), continue = FALSE, traj_death) %>% 
  ## 5.2) Did the cancer recur?
  branch(option = function() (get_attr(keys = "AdjuvantCycleEvent") == 2), continue = FALSE, traj_recurrence) %>% 
  ## 5.3) Can more treatment cycles be provided?
  rollback(amount = 5, check = function() (get_attr(keys = "AdjuvantCycles") < n_max_adjuvant_cycles)) %>% 
  
  
  ## DISEASE MONITORING AFTER ADJUVANT TREATMENT ## 
  
  # 6) Record the event and duration for the monitoring cycle
  set_attr(keys   = c("MonitorCycleEvent", "MonitorCycleTime"),
           values = function() fn_monitor_cycle(attrs = get_attr(keys = c("MonitorCycles", "BS", "TTR")),
                                                t_now = now_())) %>% 
  
  # 7) Update the number of monitoring cycles, costs, and QALYs
  set_attr(keys   = c("MonitorCycles", "dCosts", "dQALYs"),
           values = function() fn_monitor_impact(attrs = get_attr(keys = c("MonitorCycleTime")),
                                                 t_now = now_()),
           mod    = "+") %>% 
  
  # 8) Delay for the duration of the monitoring cycle
  timeout_from_attribute(key = "MonitorCycleTime") %>% 
  
  # 9) Check what will happen next based on the previously set MonitorCycleEvent:
  #     1 = death during cycle
  #     2 = recurrence during cycle
  #     3 = no recurrence or death during cycle
  ## 9.1) Did the individual die?
  branch(option = function() (get_attr(keys = "MonitorCycleEvent") == 1), continue = FALSE, traj_death) %>% 
  ## 9.2) Did the cancer recur?
  branch(option = function() (get_attr(keys = "MonitorCycleEvent") == 2), continue = FALSE, traj_recurrence) %>% 
  ## 9.3) Can more treatment cycles be provided?
  rollback(amount = 5, check = function() (get_attr(keys = "MonitorCycles") < n_max_monitor_cycles)) %>% 
  
  
  ## LONG-TERM FOLLOW UP ##
  
  # 10) Record the time until recurrence or non-cancer death
  # Here the events are defined as follows:
  #   1 = death
  #   2 = recurrence
  set_attr(keys   = c("FollowUpTime", "FollowUpEvent"),
           values = function() c(FollowUpTime  = min(get_attr(keys = c("BS", "TTR"))) - now_(),
                                 FollowUpEvent = which.min(get_attr(keys = c("BS", "TTR"))))) %>% 
  
  # 11) Update the QALYs
  set_attr(keys   = "dQALYs",
           mod    = "+",
           values = function() fn_discount_QALYs(utility    = u_diseasefree,
                                                 t_start    = now_(),
                                                 t_duration = get_attr(keys = "FollowUpTime"))) %>% 
  
  # 12) Delay until recurrence or non-cancer death
  timeout_from_attribute(key = "FollowUpTime") %>% 
  
  # 13) Was the event death or cancer recurrence?
  # There are two options here and, hence, two different trajectories the individual
  # can go to in the branch. What happens has previously been recorded in FollowUpEvent:
  #   1 = death
  #   2 = recurrence 
  branch(option = function() get_attr(keys = "FollowUpEvent"), continue = FALSE,
         
         # option = 1 -> death/BS
         traj_death,
         
         # option = 2 -> recurrence/TTR
         traj_recurrence);

# This was the end of the main trajectory.

# You can visualize the trajectory to check, for example, whether the rollback amount
# was set correctly.
plot(traj_main);




### 5. RUN SIMULATION ----

# Setting the number of individuals to be simulated:
# - 10k individuals takes about 24 seconds for both strategies
# - 100k individuals takes about 4 minutes minutes
# - 50k is sufficient for stable outcomes (see Section 6), which takes about 2 minutes
n_individuals <- 1000;

# Define the simulation environment
# - Note that the name needs to be "sim" because otherwise the short-cut functions
#   don't work properly. Other names can be used but than that name should be
#   specified in all the function calls in the trajectories or the default values
#   should be updated in the short-cut functions.
sim <- simmer() %>% 
  add_generator(name_prefix = "Patient_",
                trajectory  = traj_main,
                distribution = at(rep(x = 0, times = n_individuals)),
                mon = 2);

# Track how long the simulation takes
system.time({

  # For each treatment strategy:
  # - Set the treatment strategy
  # - Set the random number seed
  # - Reset and run the model
  # - Extract and summarize outcomes
  random_seed <- 1;
  
  # Run the simulation and summarize outcomes for strategy: Lev (1)
  treatment_arm <- 1; 
  set.seed(random_seed); 
  sim %>% reset() %>% run(); 
  df_1 <- fn_summarize(df_sim = get_mon_attributes(sim)); 

  # Run the simulation and summarize outcomes for strategy: Lev+5FU (2)
  treatment_arm <- 2; 
  set.seed(random_seed); 
  sim %>% reset() %>% run();
  df_2 <- fn_summarize(df_sim = get_mon_attributes(sim));
  
});


# Observe the mean outcomes (based on 100k individuals and random_seed = 1)
# - syntax for data.table package
# 1)  Costs     QALYs
#     76593.79  6.527653
# 2)  Costs     QALYs
#     170256.8  8.423909
df_1[ , .(Costs = mean(dCosts), QALYs = mean(dQALYs))];
df_2[ , .(Costs = mean(dCosts), QALYs = mean(dQALYs))];

# Incremental mean outcomes for 2 vs. 1
# Costs:  93663.03
# QALYS:  1.896256
mean(df_2$dCosts) - mean(df_1$dCosts);
mean(df_2$dQALYs) - mean(df_1$dQALYs);

# Incremental cost-effectiveness ratio: 49393.67
(mean(df_2$dCosts) - mean(df_1$dCosts)) / (mean(df_2$dQALYs) - mean(df_1$dQALYs));




### 6. STABILITY OF OUTCOMES ----

# Generate matrix to assess the stability of the outcomes by iteratively increasing the
# number of simulated individuals included - to generate the figure in the manuscript, 
# the simulation based on 100k individuals was used.
# - to calculate the net health benefit (NHB) a willingness to pay (WTP) of 50k is used
m_stability <- sapply(seq(1000, nrow(df_1), 1000), function(n) {
  
  # Determine the mean costs based on n (i.e., the number of individuals considered)
  costs_1 <- mean(df_1$dCosts[1:n]);
  costs_2 <- mean(df_2$dCosts[1:n]);
  
  # Determine the mean QALYs based on n (i.e., the number of individuals considered)
  QALYs_1 <- mean(df_1$dQALYs[1:n]);
  QALYs_2 <- mean(df_2$dQALYs[1:n]);
  
  # Return the incremental outcomes in a vector
  c(n     = n, 
    Costs = (costs_2 - costs_1),
    QALYs = (QALYs_2 - QALYs_1),
    ICER  = (costs_2 - costs_1) / (QALYs_2 - QALYs_1),
    NHB   = (QALYs_2 - QALYs_1) - (costs_2 - costs_1) / 50000);
  
});


# Plot the outcomes as a function of the number of individuals simulated
par(mfrow = c(2,2)); # plot in a 2x2 layout

plot(x = m_stability["n", ]/10^3, y = m_stability["Costs", ]/10^3, type = "l", las = 1,
     ylab = "Costs (thousands)", xlab = "Number of simulated individuals (thousands)",
     main = "a) Incremental Costs");

plot(x = m_stability["n", ]/10^3, y = m_stability["QALYs", ], type = "l", las = 1,
     ylab = "QALYs", xlab = "Number of simulated individuals (thousands)",
     main = "b) Incremental QALYs");

plot(x = m_stability["n", ]/10^3, y = m_stability["NHB", ], type = "l", las = 1,
     ylab = "NHB", xlab = "Number of simulated individuals (thousands)",
     main = "c) NHB (WTP = 50k)");

plot(x = m_stability["n", ]/10^3, y = m_stability["ICER", ]/10^3, type = "l", las = 1,
     ylab = "ICER (thousands)", xlab = "Number of simulated individuals (thousands)",
     main = "d) ICER");

par(mfrow = c(1,1)); # reset plot layout
