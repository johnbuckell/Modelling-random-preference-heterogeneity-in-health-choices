#####################################################################################################################
# Code written by John Buckell and Thomas Hancock (Health Economics Research Centre                                 #   
# and University of Oxford; Choice Modelling Centre) for the analysis of simulated choice                           #
# data. The simulated is based on the example drug choice data available in the Apollo package (link below).        #
#                                                                                                                   #
# The code and data are intended to help members of the community understand the choice models that we used.        #
# see also http://www.apollochoicemodelling.com/index.html                                                          #
#####################################################################################################################


# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Install Apollo (if needed)
#install.packages(apollo)

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(readr)
library(tidyverse)
library(tictoc)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName ="HIV_prev_n_u_t_ln_lu_AT_fm2_fm3_MA",# Make sure to use a new name for every model
  indivID   ="idd",  # Name of column in the database with each individual's ID
  modelDescr ="HIV_prev_n_u_t_ln_lu_AT_fm2_fm3_MA",
  mixing    = FALSE, # TRUE for models that include random parameters
  nCores    = 4,    # Number of cores to use in estimation
  seed_draws= 13    # seed for random draw generation
)

### Load data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
database = read.csv("Sim_Drug_Choice_2.csv")

database=read.csv("SimDrugChoice2_MMNL_NLN_LL.csv")
#database<-left_join(database,LL_NLN,by="id")
#rm(LL_NLN)
LL_U=read.csv("SimDrugChoice2_MMNL_uniforms_LL.csv")
database<-left_join(database,LL_U,by="idd")
rm(LL_U)
LL_T=read.csv("SimDrugChoice2_MMNL_triangulars_LL.csv")
database<-left_join(database,LL_T,by="idd")
rm(LL_T)
LL_LN=read.csv("SimDrugChoice2_MMNL_lognormal_LL.csv")
database<-left_join(database,LL_LN,by="idd")
rm(LL_LN)
LL_LU=read.csv("SimDrugChoice2_MMNL_loguniform_LL.csv")
database<-left_join(database,LL_LU,by="idd")
rm(LL_LU)
LL_AT=read.csv("SimDrugChoice2_AT_LL.csv")
database<-left_join(database,LL_AT,by="idd")
rm(LL_AT)
LL_fm2=read.csv("SimDrugChoice2_MMNL_FM2_LL.csv")
database<-left_join(database,LL_fm2,by="idd")
rm(LL_fm2)
LL_fm3=read.csv("SimDrugChoice2_MMNL_FM3_LL.csv")
database<-left_join(database,LL_fm3,by="idd")
rm(LL_fm3)
database<-as.data.frame(database)


# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(

  delta_normal =0,
  delta_uniform =3.608,
  delta_triangular =-2.9,
  delta_lognormal =-8.669,
  delta_loguniform =3.817,
  delta_at =-2.937,
  delta_fm2 =5.764,
  delta_fm3 =7.253
)
apollo_fixed = c("delta_normal")

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none

apollo_lcPars=function(apollo_beta, apollo_inputs){
  lcpars = list()
  
  #### Model averaging MNL
  Vdr=list()
  Vdr[["class_normal"]] = delta_normal
  Vdr[["class_uniform"]] = delta_uniform
  Vdr[["class_triangular"]] = delta_triangular
  Vdr[["class_lognormal"]] = delta_lognormal
  Vdr[["class_loguniform"]] = delta_loguniform
  Vdr[["class_at"]] = delta_at
  Vdr[["class_fm2"]] = delta_fm2
  Vdr[["class_fm3"]] = delta_fm3
  eVdr <<- lapply(Vdr,exp)
  

  # Calculate class membership probabilities using MNL specification
  lcpars[["pi_values"]][["class_normal"]]=(eVdr[["class_normal"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_uniform"]]=(eVdr[["class_uniform"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_triangular"]]=(eVdr[["class_triangular"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_lognormal"]]=(eVdr[["class_lognormal"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_loguniform"]]=(eVdr[["class_loguniform"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_at"]]=(eVdr[["class_at"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_fm2"]]=(eVdr[["class_fm2"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]][["class_fm3"]]=(eVdr[["class_fm3"]]/(eVdr[["class_normal"]]+eVdr[["class_uniform"]]+eVdr[["class_triangular"]]+eVdr[["class_lognormal"]]+eVdr[["class_loguniform"]]+eVdr[["class_at"]]+eVdr[["class_fm2"]]+eVdr[["class_fm3"]]))
  lcpars[["pi_values"]] = apollo_firstRow(lcpars[["pi_values"]], apollo_inputs)
  test<<-lcpars
  
  return(lcpars)
}


apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ## normal
  P[['mnl_n']] = exp(avgLL_NLN)
  
  ## uniform
  P[['mnl_u']] = exp(avgLL_U)
  
  ## triangular
  P[['mnl_t']] = exp(avgLL_T)

  ## lognormal
  P[['mnl_ln']] = exp(avgLL_LN)
  
  ## loguniform
  P[['mnl_lu']] = exp(avgLL_LU)
  
  ## at
  P[['mnl_at']] = exp(avgLL_AT)
  
  ## fm2
  P[['mnl_fm2']] = exp(avgLL_FM2)
  
  ## fm3
  P[['mnl_fm3']] = exp(avgLL_FM3)
  
  ### Compute latent class model probabilities
  lc_settings   = list(inClassProb = P, classProb=pi_values)
  P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
  
}

# ################################################################# #
##### ESTIMATE THE MODEL AND GET OUTPUT                          ####
# ################################################################# #
settings<-list(constraints = c(),estimationRoutine="bgw",maxIterations=500)

### Estimate model
model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,estimate_settings=settings)
apollo_modelOutput(model)

