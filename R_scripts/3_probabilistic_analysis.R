# This script belongs to the manuscript: Degeling K, Karnon J, Van de Ven M, Brennan A,
# Koffijberg H. Discrete event simulation in R using the 'simmer' package for health
# economic modeling: a tutorial and illustration in colon cancer. Applied Health 
# Economics and Health Policy.
# 
# This is the script that performs the probabilistic analysis of the discrete event
# simulation (DES) defined in the 2_deterministic_model.R script.
#
# The script has the following sections:
# 1. INITIALIZATION               loading the packages and survival analysis results
# 2. SAMPLING PARAMETER VALUES    sampling parameter values for use in the analysis
# 3. SIMULATION FUNCTION          wrapping the DES in a function
# 4. PROBABILISTIC ANALYSIS       performing the probabilistic analysis of the DES
# 5. OUTPUT ANALYSIS              analyzing the output of the probabilistic analysis
#
# Version history:
# v1    28-10-2020  Koen Degeling     R 3.6.1   Initial version
# v2    20-07-2021  Koen Degeling     R 4.0.3   Updating to use FORK rather than SOCK 
#                                               cluster, and use with() function rather
#                                               than loading parameters separately




### 1. INITIALIZATION ----

# Uncomment to clear objects from Global Environment and clear the console
#rm(list = ls()); cat("\014");

# Required packages
library(simmer);        # version 4.4.2
library(simmer.plot);   # version 0.1.16
library(flexsurv);      # version 2.0 - for sampling from distributions
library(data.table);    # version 1.14.0 - for summarizing the simulation output
library(parallel);      # version 4.0.3 - for running the analysis on multiple cores

# Loading the parameters estimated in the survival analysis
# - only used for background survival (BS)
df_survival_models <- read.csv(file = "data/survival_models.csv");

# Loading the model fits from the survival analysis
# - used to reflect uncertainty in the parameters for time-to-recurrence (TTR) and time
#   from recurrence to death due to cancer (TTD)
load(file = "data/survival_models.RDS");




### 2. SAMPLING PARAMETER VALUES ----

# This section generates the parameter value samples that will be used to perform the
# probabilistic analysis. The values used to define the parametric distributions that
# are used to generate the samples are based on assumptions as highlighted in the 
# corresponding comments. Only for time-to-recurrence and time-to-death from recurrence,
# the parameter uncertainty is based on the data used in the survival analysis.

# Set the number of random samples that are to be generated for each parameter and set
# a random number seed for reproducibility
n_samples <- 1000;
set.seed(1);

# Initialise a data.frame that will store all the parameter samples
df_pa <- data.frame(i_run = 1:n_samples);

# Probability of being high risk.
# - no uncertainty considered, because this was selected as population of interest
df_pa$p_highrisk <- 0.6;

# Gompertz distribution for death due to other causes
# - no uncertainty considered, because modeled based on life tables
df_pa$d_BS_shape <- df_survival_models$d_RF_D_shape;
df_pa$d_BS_rate  <- df_survival_models$d_RF_D_rate;

# Log-logistic cure distribution for time-to-recurrence (TTR)
# - correlated sets of parameter values are generated using multivariate normal 
#   distributions using the mvrnorm() function of the base R MASS package.
# - note that these parameters are not always on the right scale, so some
#   transformations are required - see 1_survival_analysis.R for details.
m_TTR <- MASS::mvrnorm(n = n_samples, mu = fit_RF_R$coefficients, Sigma = fit_RF_R$cov);

thetaToProb <- function(x) unname(exp(x) / (1 + exp(x)));
df_pa$d_TTR_cure_1_low  <- thetaToProb(m_TTR[, "theta"]);
df_pa$d_TTR_cure_1_high <- thetaToProb(m_TTR[, "theta"] + m_TTR[, "node4"]);
df_pa$d_TTR_cure_2_low  <- thetaToProb(m_TTR[, "theta"] + m_TTR[, "rxLev+5FU"]);
df_pa$d_TTR_cure_2_high <- thetaToProb(m_TTR[, "theta"] + m_TTR[, "rxLev+5FU"] + m_TTR[, "node4"]);

df_pa$d_TTR_shape_1_low  <- exp(m_TTR[, "shape"]);
df_pa$d_TTR_shape_1_high <- exp(m_TTR[, "shape"]);
df_pa$d_TTR_shape_2_low  <- exp(m_TTR[, "shape"]);
df_pa$d_TTR_shape_2_high <- exp(m_TTR[, "shape"]);

df_pa$d_TTR_scale_1_low  <- exp(m_TTR[, "scale"]);
df_pa$d_TTR_scale_1_high <- exp(m_TTR[, "scale"] + m_TTR[, "scale(node4)"]);
df_pa$d_TTR_scale_2_low  <- exp(m_TTR[, "scale"] + m_TTR[, "scale(rxLev+5FU)"]);
df_pa$d_TTR_scale_2_high <- exp(m_TTR[, "scale"] + m_TTR[, "scale(rxLev+5FU)"] + m_TTR[, "scale(node4)"]);

# Adjuvant treatment costs
# - no uncertainty considered, because assumed to be fixed
df_pa$c_adjuvant_cycle_1 <- 5000;
df_pa$c_adjuvant_cycle_2 <- 15000;

# Utility during adjuvant treatment
# - assume: mean = 0.7, se = 0.7*0.001
df_pa$u_adjuvant <- rlnorm(n = n_samples, meanlog = -0.35738872, sdlog = 0.03778296);

# Probability of toxicities during a cycle
# - assume: based on 500 individuals in each arm and p_1 = 0.1 and p_2 = 0.2
df_pa$p_tox_1 <- rbeta(n = n_samples, shape1 = 0.1*500, shape2 = (1-0.1)*500);
df_pa$p_tox_2 <- rbeta(n = n_samples, shape1 = 0.2*500, shape2 = (1-0.2)*500);

# Costs of toxicities
# - assume: mean = 2000, se = 2000*0.2
df_pa$c_tox <- rgamma(n = n_samples, shape = 10000, rate = 5);

# Disutility of toxicities
# - assume: mean = 0.1, se = 0.1*0.001
df_pa$u_dis_tox <- rlnorm(n = n_samples, meanlog = -2.30756026, sdlog = 0.09975135);

# Other adjuvant treatment parameters
# - no uncertainty
df_pa$t_adjuvant_cycle      <- 3/52;
df_pa$n_max_adjuvant_cycles <- 10;

# Cost of a monitoring cycle
# - no uncertainty considered, because assumed to be fixed
df_pa$c_monitor <- 1000;

# Utility while free of recurrence
# - assume: mean = 0.8, se = 0.8*0.001
df_pa$u_diseasefree <- rlnorm(n = n_samples, meanlog = -0.2237682, sdlog = 0.0353443);

# Other disease monitoring parameters
df_pa$t_monitor_cycle      <- 1; # in years
df_pa$n_max_monitor_cycles <- 5;

# Log-logistic distribution for time from recurrence to death due to cancer (TTD)
# - correlated sets of parameter values are generated using multivariate normal 
#   distributions using the mvrnorm() function of the base R MASS package.
# - note that these parameters are not always on the right scale, so some
#   transformations are required - see 1_survival_analysis.R for details.
m_TTD <- MASS::mvrnorm(n = n_samples, mu = fit_R_D$coefficients, Sigma = fit_R_D$cov);

df_pa$d_TTD_shape_1 <- exp(m_TTD[, "shape"]);
df_pa$d_TTD_shape_2 <- exp(m_TTD[, "shape"] + m_TTD[, "shape(rxLev+5FU)"]);

df_pa$d_TTD_scale_1 <- exp(m_TTD[, "scale"]);
df_pa$d_TTD_scale_2 <- exp(m_TTD[, "scale"] + m_TTD[, "rxLev+5FU"]);

# Cost of treatment for advanced disease
# - assume: mean = 40000, se = 40000*0.2
df_pa$c_advanced <- rgamma(n = n_samples, shape = 200000, rate = 5);

# Utility during advanced disease
# - assume: mean = 0.6, se = 0.6*0.001
df_pa$u_advanced <- rlnorm(n = n_samples, meanlog = -0.51165826, sdlog = 0.04080783);




### 3. SIMULATION FUNCTION ----

# The runDES() function encapsulates basically the complete 1_deterministic_model.R
# script and, thereby, can be used to perform a run of the model.
runDES <- function(parameters, n_individuals = 1000, summarize_output = TRUE, random_seed = 1, dr_costs = 0.04, dr_health = 0.015) {
  
  # Using the with() function, the code within the expression (expr = {...}) can be run
  # without the need to extract all the parameter values from the parameters object,
  # because the names in "parameters" are in line with those in the text. In other
  # words, when a parameter is used in the code, R will look in the "parameters" object
  # for a parameter with that name and use that.
  with(data = parameters, expr = {
    
    # Store the time at which the simulation is started to track the run time.
    tic <- Sys.time();
    
    
    ## 2.1 SUPPORTING FUNCTIONS ----
    
    # The supporting functions need to be defined in the same environment as the simmer
    # trajectory, otherwise it does not know where to find them. Alternatively to
    # defining them below, they can be placed in a separate script and "sourced" (i.e.,
    # loaded) here using the source() function.
    # - see the 2_deterministic_model.R script for details about the functions
    
    fn_discount_costs <- function(amount, t) amount / (1 + dr_costs) ^ t;
    fn_discount_QALYs <- function(utility, t_start, t_duration) utility / (1 + dr_health) ^ t_start / log(1 + dr_health) * (1 - 1 / (1 + dr_health) ^ t_duration);
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
    rllogiscure  <- function(n = 1, cure, shape, scale) ifelse(runif(n) < cure, Inf, rllogis(n = n, shape = shape, scale = scale));
    set_attr <- set_attribute;
    get_attr <- function(keys, .env = sim) get_attribute(.env = .env, keys = keys);
    now_     <- function(.env = sim) now(.env = .env);
    fn_initialisation  <- function() {
      
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
    fn_adjuvant_cycle  <- function(attrs, t_now) {
      
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
    fn_monitor_cycle   <- function(attrs, t_now) {
      
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
    fn_monitor_impact  <- function(attrs, t_now) {
      
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
    fn_advanced_time   <- function(attrs) {
      
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
    
    
    ## 2.2 MODEL STRUCTURE ----
    
    # Similar to the supporting functions, the trajectory can also be saved in a
    # separate script and sourced here using the source() function.
    
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
    
    
    ## 2.3 RUN SIMULATION ----
    
    # Define the environment
    # - note that n_individuals and random_seed have been provided as an argument to 
    #   the runDES() function
    sim <- simmer() %>% 
      add_generator(name_prefix = "Patient_",
                    trajectory  = traj_main,
                    distribution = at(rep(x = 0, times = n_individuals)),
                    mon = 2);
    
    # Simulation for strategy 1: Lev
    treatment_arm <- 1; set.seed(random_seed); sim %>% reset() %>% run(); 
    df_sim_1 <- get_mon_attributes(sim);
    
    # Simulation for strategy 2: Lev+5FU
    treatment_arm <- 2; set.seed(random_seed); sim %>% reset() %>% run();
    df_sim_2 <- get_mon_attributes(sim);
    
    
    ## 2.4 EXTRACT AND RETURN OUTCOMES ----
    
    # Conditional on the value of the summarize_output argument to the runDES() function,
    # either the complete simulation data.frame is returned (summarize_output == FALSE)
    # or the summarized health economic outcomes are returned
    if(summarize_output) {
      
      # Summarize attributes - only extract those required for efficiency
      df_1 <- fn_summarize(df_sim = df_sim_1, keys = c('AdjuvantCycles', 'Toxicities', 'MonitorCycles', 'TTR', 'TTD', 'OS', 'dCosts', 'dQALYs')); 
      df_2 <- fn_summarize(df_sim = df_sim_2, keys = c('AdjuvantCycles', 'Toxicities', 'MonitorCycles', 'TTR', 'TTD', 'OS', 'dCosts', 'dQALYs')); 
      
      # Summarize outcomes
      output <- c(
        unlist(parameters),
        AdjuvantCycles_1 = mean(df_1$AdjuvantCycles),
        AdjuvantCycles_2 = mean(df_2$AdjuvantCycles),
        Toxicities_1     = mean(df_1$Toxicities),
        Toxicities_2     = mean(df_2$Toxicities),
        MonitorCycles_1  = mean(df_1$MonitorCycles),
        MonitorCycles_2  = mean(df_2$MonitorCycles),
        Recurrence_1     = mean(!is.na(df_1$TTD)),
        Recurrence_2     = mean(!is.na(df_2$TTD)),
        TTR_1    = mean(df_1$TTR[!is.na(df_1$TTD)]),
        TTR_2    = mean(df_2$TTR[!is.na(df_2$TTD)]),
        OS_1     = mean(df_1$OS),
        OS_2     = mean(df_2$OS),
        dQALYs_1 = mean(df_1$dQALYs),
        dQALYs_2 = mean(df_2$dQALYs),
        dCosts_1 = mean(df_1$dCosts),
        dCosts_2 = mean(df_2$dCosts)
      );
      
    } else {
      
      # Summarize all attributes
      df_1 <- fn_summarize(df_sim = df_sim_1); 
      df_2 <- fn_summarize(df_sim = df_sim_2); 
      
      # Combine in single data.frame
      output <- rbind(df_1, df_2);
      
    }
    
    # Determine and print the run time
    toc <- Sys.time();
    print(toc - tic);
    
    # IMPORTANT! clean up the environment so remaining objects don't consume memory
    rm(list = ls()[ls() != "output"]); gc();
  
    return(output);
    
  }) # end of: with(data = parameters, expr = {..})
  
}

# Using the runDES function, a single run of the model for a certain set of parameter
# values can be performed by uncommenting the line below. For example, you can change
# the number of individuals to be simulated (illustrated) or another argument.
# This may be useful for testing purposes.
#runDES(parameters = df_pa[1, ]);
#runDES(parameters = df_pa[1, ], n_individuals = 10000);




### 4. PROBABILISTIC ANALYSIS ----

# IMPORTANT! The code below runs the probabilistic analysis in parallel to decrease the
# run time. This means that if your PC has 10 cores and all these are used, 10 instances
# of the model run at the same time, which may require substantial memory (up to 1.5GB
# per instance for 50 thousand simulated individuals). Therefore, it is important to
# start with a low number of individuals (and potentially cores) to see how much your
# PC can handle. To generate the results for the manuscript, a MacBook Pro 2020 with
# 12 cores (actually 6 cores, but 12 threads) and 16GB of memory, which could easily
# handle a parallel simulation using 12 threads for 50k individuals per run.
# - you can check the utilization of your PC memory using the Task Manager (Windows) or
#   the Activity Monitor (MacOS).
# - results of a probabilistic analysis of 1000 runs of 50000 individuals per run, are
#   saved in the probabilistic_analysis.RData, which can be loaded optional in Section
#   5 so that you don't have to run the analysis yourself to generate the plots.

# Specify the number of runs to be performed. Note that this number cannot be higher
# than the number of rows in the df_pa data.frame, because a set of parameters is
# required for every run.
n_runs <- 10;

# Specify the number of cores that are not to be used
# - you can use detectCores() to find out how many are available
#detectCores();
n_free_cores <- 2;

# Note that there are different ways that computing cluster can work. For the purpose of
# running a probabilistic analysis the main difference is the memory management. In a
# "fork" cluster, each thread is run separately but they use the same environment. In a
# "socket" cluster, each thread has its own environment with all objects required to 
# execute the code, which will have to be exported from the main global environment. 
# Therefore, a fork cluster is easier and quicker to set up and more efficient in terms
# of memory usage, and generally preferable. However, fork clusters do not work on
# Windows machines. Hence, below we show the implementation for both fork and socket
# clusters.


## 4.1 Fork cluster (MacOS) ----

# Define the FORK computing cluster
cl_fork <- makeForkCluster(detectCores() - n_free_cores);

# The probabilistic analysis is run in parallel using the parSapply() function, which
# is the parallel version of the sapply() function. It nicely formats the output in a
# matrix with a column for every run. In each run, a different set of parameters is used
# by forwarding a different row of the df_pa data.frame into the runDES() function. Also
# note that the default number of individuals can be overwritten here.
system.time({
  m_output <- parSapply(cl = cl_fork, X = 1:n_runs, FUN = function(i_run) runDES(parameters = df_pa[i_run, ], n_individuals = 50000));
})

# After the simulation is finished, nicely disable the computing cluster
stopCluster(cl = cl_fork);

# To analyze and plot the outcomes in the next section, the matrix is transformed to a 
# data.frame with the outcomes as columns and the runs as rows.
df_output <- as.data.frame(t(m_output));

# Uncomment this line to overwrite the saved probabilistic analysis outcomes
#save(df_output, file = "data/probabilistic_analysis.RData");


## 4.2 Socket cluster (Windows/MacOS) ----

# Define the SOCKET computing cluster
cl_sock <- makePSOCKcluster(detectCores() - n_free_cores);

# Because this is a SOCKET cluster, we will need to export all the object form the
# global environment that should be available to the threads to execute the code
# - Note that this also includes functions of specific packages, as each thread basically
#   is a new clean instance of R. To prevent exporting all the functions from, for
#   example the simmer package, one by one, we can load the required packages before
#   running the simulation as demonstrated below. That means we only need to export the
#   objects and functions that are not included in those packages.
clusterExport(cl = cl_sock, varlist = c("df_pa", "runDES"));

# As explained above, for the SOCKET cluster, we need to load the packages before we can 
# run the runDES() function, otherwise errors will be returned that certain functions
# could not be found.
system.time({
  m_output <- parSapply(cl = cl_sock, X = 1:n_runs, FUN = function(i_run) {
    
    library(simmer);
    library(flexsurv);
    library(data.table);
    
    runDES(parameters = df_pa[i_run, ], n_individuals = 50000);
  
  })
})

# After the simulation is finished, nicely disable the computing cluster
stopCluster(cl = cl_sock);

# To analyze and plot the outcomes in the next section, the matrix is transformed to a 
# data.frame with the outcomes as columns and the runs as rows.
df_output <- as.data.frame(t(m_output));

# Uncomment this line to overwrite the saved probabilistic analysis outcomes
#save(df_output, file = "data/probabilistic_analysis.RData");




### 5. OUTPUT ANALYSIS ----

# Uncomment the line below to load the saved probabilistic analysis results.
#load(file = "data/probabilistic_analysis.RData");

# Generate matrix to assess the stability of the outcomes by iteratively increasing the
# number of probabilistic runs included
# - to calculate the net health benefit (NHB) a willingness to pay (WTP) of 50k is used
m_stability <- sapply(X = seq(from = 10, to = 1000, by = 10), FUN = function(n) {
  
  # Extract the runs to be included
  df <- df_output[1:n, ];
  
  # Incremental costs including 95% confidence interval
  inc_costs_mean <- mean(df$dCosts_2 - df$dCosts_1);
  inc_costs_lb   <- unname(quantile(df$dCosts_2 - df$dCosts_1, 0.025));
  inc_costs_ub   <- unname(quantile(df$dCosts_2 - df$dCosts_1, 0.975));
  
  # Incremental QALYs including 95% confidence interval
  inc_QALYs_mean <- mean(df$dQALYs_2 - df$dQALYs_1);
  inc_QALYs_lb   <- unname(quantile(df$dQALYs_2 - df$dQALYs_1, 0.025));
  inc_QALYs_ub   <- unname(quantile(df$dQALYs_2 - df$dQALYs_1, 0.975));  
  
  # Net health benefit including 95% confidence interval
  NHB_mean       <- mean((df$dQALYs_2 - df$dQALYs_1) - (df$dCosts_2 - df$dCosts_1)/50000);
  NHB_lb         <- unname(quantile((df$dQALYs_2 - df$dQALYs_1) - (df$dCosts_2 - df$dCosts_1)/50000, 0.025));
  NHB_ub         <- unname(quantile((df$dQALYs_2 - df$dQALYs_1) - (df$dCosts_2 - df$dCosts_1)/50000, 0.975));
  
  # Return outcomes
  c(n = n,
    inc_costs_mean = inc_costs_mean,
    inc_costs_lb   = inc_costs_lb,
    inc_costs_ub   = inc_costs_ub,
    inc_QALYs_mean = inc_QALYs_mean,
    inc_QALYs_lb   = inc_QALYs_lb,
    inc_QALYs_ub   = inc_QALYs_ub,
    NHB_mean = NHB_mean,
    NHB_lb   = NHB_lb,
    NHB_ub   = NHB_ub
  );
  
});

# Plot the outcomes as a function of the number of probabilistic runs
# - 1000 is more than enough
plot(x = m_stability["n", ], y = m_stability["inc_costs_mean", ]/1000,
     ylim = c(80, 100), type = "l", las = 1, 
     main = "Incremental Costs",
     xlab = "Number of probabilistic analysis runs", 
     ylab = "Costs (thousands)");
lines(x = m_stability["n", ], y = m_stability["inc_costs_lb", ]/1000, lty = 2);
lines(x = m_stability["n", ], y = m_stability["inc_costs_ub", ]/1000, lty = 2);

plot(x = m_stability["n", ], y = m_stability["inc_QALYs_mean", ],
     ylim = c(0, 4), type = "l", las = 1, 
     main = "Incremental QALYs",
     xlab = "Number of probabilistic analysis runs", 
     ylab = "QALYs");
lines(x = m_stability["n", ], y = m_stability["inc_QALYs_lb", ], lty = 2);
lines(x = m_stability["n", ], y = m_stability["inc_QALYs_ub", ], lty = 2);

plot(x = m_stability["n", ], y = m_stability["NHB_mean", ],
     ylim = c(-2, 2), type = "l", las = 1, 
     main = "NHB (WTP = 50k/QALY)",
     xlab = "Number of probabilistic analysis runs", 
     ylab = "NHB");
lines(x = m_stability["n", ], y = m_stability["NHB_lb", ], lty = 2);
lines(x = m_stability["n", ], y = m_stability["NHB_ub", ], lty = 2);


# Observe the mean outcomes and confidence intervals:
# 1)  Costs                           QALYs
#     74685.27 (72586.43, 76650.75)   6.557272 (5.816502, 7.367419)
# 2)  Costs                           QALYs
#     166438.21 (163340.4, 169249.3)  8.438214 (7.463922, 9.475895)
colMeans(df_output[, c("dCosts_1", "dCosts_2")]);
colMeans(df_output[, c("dQALYs_1", "dQALYs_2")]);

quantile(df_output[, "dCosts_1"], probs = c(0.025, 0.975));
quantile(df_output[, "dCosts_2"], probs = c(0.025, 0.975));

quantile(df_output[, "dQALYs_1"], probs = c(0.025, 0.975));
quantile(df_output[, "dQALYs_2"], probs = c(0.025, 0.975));

# Incremental mean outcomes for 2 vs. 1 and confidence intervals
# Costs:  91752.94 (88187.06, 95105.39)
# QALYS:  1.880942 (0.8872122, 2.9545004)
mean(df_output$dCosts_2) - mean(df_output$dCosts_1);
mean(df_output$dQALYs_2) - mean(df_output$dQALYs_1);

quantile((df_output$dCosts_2 - df_output$dCosts_1), probs = c(0.025, 0.975));
quantile((df_output$dQALYs_2 - df_output$dQALYs_1), probs = c(0.025, 0.975));

# Incremental cost-effectiveness ratio for 2 vs. 1: 48780.32
(mean(df_output$dCosts_2) - mean(df_output$dCosts_1)) / (mean(df_output$dQALYs_2) - mean(df_output$dQALYs_1));

# Net health benefit for 2 vs. 1 and confidence interval (WTP = 50k/QALY): 0.04588299 (-0.9870472, 1.1379313)
(mean(df_output$dQALYs_2) - mean(df_output$dQALYs_1)) - (mean(df_output$dCosts_2) - mean(df_output$dCosts_1))/50000;

quantile((df_output$dQALYs_2 - df_output$dQALYs_1) - (df_output$dCosts_2 - df_output$dCosts_1)/50000, probs = c(0.025, 0.975));


# Incremental cost-effectiveness plane
plot(x = c(0, 10), y = c(0, 10*50),
     xlim = c(0, 4), ylim = c(70, 100),
     type = "l", col = "red", lwd = 2, las = 1,
     xlab = "Incremental QALYs", ylab = "Incremental Costs (thousands)",
     main = "Incremental Cost-Effectiveness Plane");
car::dataEllipse(x = (df_output$dQALYs_2 - df_output$dQALYs_1),
                 y = (df_output$dCosts_2 - df_output$dCosts_1)/1000,
                 grid = FALSE, levels = 0.95, add = TRUE);
legend("bottomright", bty = "n",
       legend = c("Probabilistic Runs", "Mean Outcomes", "95% Confidence Ellipse", "WTP (50k/QALY)"),
       col = c("black", "blue", "blue", "red"),
       lty = c(NA, NA, 1, 1),
       pch = c(1, 16, NA, NA));


# Create a summary table of the probabilistic analysis:
#   1) Add a column with the differences
#   2) Select the columns to be included (i.e., excluding the parameter values)
#   3) Pivot to long format
#   4) Extract the strategy (or whether it is a difference) from the variable name
#   5) Calculate the summary statistics (mean and 95% confidence interval)
#   6) Combine the summary statistics into a easy to read format: est (lb, ub)
#   7) Pivot to wider format to obtain the table
tbl_output <- data.table(df_output)[ , `:=` (
  AdjuvantCycles_d = AdjuvantCycles_2 - AdjuvantCycles_1,
  Toxicities_d     = Toxicities_2 - Toxicities_1,
  MonitorCycles_d  = MonitorCycles_2 - MonitorCycles_1,
  Recurrence_d     = Recurrence_2 - Recurrence_1,
  TTR_d    = TTR_2 - TTR_1,
  OS_d     = OS_2 - OS_1,
  dCosts_d = dCosts_2 - dCosts_1,
  dQALYs_d = dQALYs_2 - dQALYs_1,
  NHB_1    = dQALYs_1 - dCosts_1/50000,
  NHB_2    = dQALYs_2 - dCosts_2/50000
)][ , NHB_d := NHB_2 - NHB_1]
tbl_output <- tbl_output[, .SD, .SDcols = names(tbl_output) %like% 'AdjuvantCycles|Toxicities|MonitorCycles|Recurrence|^TTR|OS|dCosts|dQALYs|NHB']
tbl_output <- melt(tbl_output, id.vars = NULL, measure.vars = names(tbl_output), variable.factor = FALSE)
tbl_output[ , `:=` (
  strategy = substr(x = variable, start = nchar(variable), stop = nchar(variable)),
  variable = substr(x = variable, start = 1, stop = nchar(variable)-2)
)]
tbl_output <- tbl_output[ , .(est = mean(value), lb = quantile(value, probs = 0.025), ub = quantile(value, probs = 0.975)), by = list(strategy, variable)]
tbl_output[ , value := paste0(round(est, 2), ' (', round(lb, 2), ', ', round(ub, 2), ')')]
tbl_output <- dcast(data = tbl_output, formula = variable~strategy)

# Save as CSV file
write.csv(x = tbl_output, file = 'data/probabilistic_analysis_output.csv', row.names = FALSE)

