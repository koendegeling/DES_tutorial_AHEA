DISCRETE EVENT SIMULATION IN R USING THE SIMMER PACKAGE FOR HEALTH ECONOMIC MODELING: A TUTORIAL AND ILLUSTRATION IN COLON CANCER

This repository includes all scripts and objects related to the tutorial on modeling health and economic outcomes using discrete event simulation in R using the simmer package. Carefully read this brief README file before running the code in the scripts.

The following important folders and files are included in the repository:

DES_tutorial_AHEA
  |
  +-- R-scripts
  |   |
  |   +-- 1_survival_analysis.R         Script including the survival analysis
  |   +-- 2_deterministic_model.R       Script for the deterministic analysis
  |   \-- 3_probabilistic_analysis.R    Script for the probabilstic analysis
  |
  |-- data
  |   |
  |   +-- probabilistic_analysis.RData  Saved result of the probabilistic analysis
  |   +-- survival_models.RDS           Survival analysis objects
  |   \-- survival_models.csv           Survival analysis parameter values
  |
  +-- DES_tutorial_AHEA.Rproj           R Studio project file

IMPORTANT! To run all scripts, ensure the working directory in R is set according to the DES_tutorial_AHEA folder as illustrated above. A different name for the main folder may be used, but it should include the data folder with the corresponding files. This can be easily be achieved by opening the DES_tutorial_AHEA.Rproj project file in R Studio.

The citation of the tutorial is as follows:
Degeling, K., Karnon, J., van de Ven, M. et al. Discrete Event Simulation in R using the ‘Simmer’ Package for Health Economic Modelling: A Tutorial and Illustration in Colon Cancer. Appl Health Econ Health Policy (2025). https://doi.org/10.1007/s40258-025-00983-8
