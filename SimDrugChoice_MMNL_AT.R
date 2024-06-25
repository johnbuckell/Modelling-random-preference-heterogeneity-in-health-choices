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

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(readr)
library(tidyverse)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "SimDrugChoice2_AT",
  modelDescr      = "SimDrugChoice2_AT",
  indivID         = "ID", 
  outputDirectory = "output",
  noDiagnostics   = TRUE,
  #noValidation    = TRUE, 
  mixing    = TRUE, # TRUE for models that include random parameters
  nCores    = 4,    # Number of cores to use in estimation
  seed_draws= 13    # seed for random draw generation
)

### Load data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
database = read.csv("Sim_Drug_Choice_2.csv")

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(
  b_unbranded =0,
  b_char_standard =0,
  branded_a =-5.795375,
  branded_b =5.10275,
  branded_c_delta =0.51,
  country_USA_bound =0,
  country_USA_range_t =1.684235,
  country_CH_bound =-0.749596,
  country_CH_range_t =2.470574,
  country_DK_bound =0.210414,
  country_DK_range_t =1.458513,
  country_IND_bound =-2.229137,
  country_IND_range_t =3.993151,
  country_RUS_bound =0.513711,
  country_RUS_range_t =1.241524,
  country_BRA_bound =-1.038217,
  country_BRA_range_t =2.697997,
  char_fast_bound =0.223935,
  char_fast_range_t =0.791833,
  char_double_bound =0.938052,
  char_double_range_t =1.065621,
  risk_bound =-0.003133,
  risk_range_t =0.001221,
  price_a =-1.479558,
  price_b =0.052078,
  price_c_delta =0.476
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed=names(apollo_beta)[apollo_beta==0|apollo_beta==1]

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "mlhs", # use mlhs when using 5 or more distributions 
  interNDraws    = 500,
  interUnifDraws = c("draws_unif1_branded1",
                     "draws_unif1_branded2",
                     "draws_unif1_CH",
                     "draws_unif1_USA",
                     "draws_unif1_RUS",
                     "draws_unif1_BRA",
                     "draws_unif1_DK",
                     "draws_unif1_IND",
                     "draws_unif1_fast",
                     "draws_unif1_double",
                     "draws_unif1_risk",
                     "draws_unif1_price",
                     
                     "draws_unif2_branded1",
                     "draws_unif2_branded2",
                     "draws_unif2_CH",
                     "draws_unif2_USA",
                     "draws_unif2_RUS",
                     "draws_unif2_BRA",
                     "draws_unif2_DK",
                     "draws_unif2_IND",
                     "draws_unif2_fast",
                     "draws_unif2_double",
                     "draws_unif2_risk",
                     "draws_unif2_price"
  ),
  interNormDraws = c("draws_norm_branded",
                     "draws_norm_CH",
                     "draws_norm_USA",
                     "draws_norm_RUS",
                     "draws_norm_BRA",
                     "draws_norm_DK",
                     "draws_norm_IND",
                     "draws_norm_fast",
                     "draws_norm_double",
                     "draws_norm_risk",
                     "draws_norm_price"
  ),
  
  intraDrawsType = "mlhs",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)



### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  randcoeff[["branded_lower"]] = branded_a + (((branded_a+branded_b)/2 + branded_c_delta) - branded_a) * sqrt(draws_unif1_branded1) 
  randcoeff[["branded_upper"]] = branded_b - (branded_b - ((branded_a+branded_b)/2 + branded_c_delta)) * sqrt(draws_unif2_branded1) 
  
  randcoeff[["price_lower"]] = price_a + (((price_a+price_b)/2 + price_c_delta) - price_a) * sqrt(draws_unif1_branded2) 
  randcoeff[["price_upper"]] = price_b - (price_b - ((price_a+price_b)/2 + price_c_delta)) * sqrt(draws_unif2_branded2) 
  
  randcoeff[["b_country_USA"]]  = country_USA_bound + country_USA_range_t * (draws_unif1_USA + draws_unif2_USA) 
  randcoeff[["b_country_CH"]]  = country_CH_bound + country_CH_range_t * (draws_unif1_CH + draws_unif2_CH) 
  randcoeff[["b_country_DK"]]  = country_DK_bound + country_DK_range_t * (draws_unif1_DK + draws_unif2_DK) 
  randcoeff[["b_country_IND"]] = country_IND_bound + country_IND_range_t * (draws_unif1_IND + draws_unif2_IND) 
  randcoeff[["b_country_RUS"]] = country_RUS_bound + country_RUS_range_t * (draws_unif1_RUS + draws_unif2_RUS) 
  randcoeff[["b_country_BRA"]] = country_BRA_bound + country_BRA_range_t * (draws_unif1_BRA + draws_unif2_BRA) 
  randcoeff[["b_char_fast"]]   = char_fast_bound + char_fast_range_t * (draws_unif1_fast + draws_unif2_fast) 
  randcoeff[["b_char_double"]] = char_double_bound + char_double_range_t * (draws_unif1_double + draws_unif2_double) 
  randcoeff[["b_risk"]]        = risk_bound + risk_range_t * (draws_unif1_risk + draws_unif2_risk) 
  
  
  return(randcoeff)
}


# ################################################################# #
#### DEFINE LATENT CLASS COMPONENTS                              ####
# ################################################################# #
lc_set <- expand.grid(
  "branded"=c("branded_lower","branded_upper"), 
  "risk"=c("price_lower","price_upper")
)

lc_set2 <- expand.grid(
  "branded"=c("branded_lower","branded_upper"), 
  "risk"=c("price_lower","price_upper")
)


apollo_lcPars=function(apollo_beta, apollo_inputs){
  lcpars = list()
  
  
  ### Create empty lists for parameters in classes and class allocation probabilities
  lcpars[["branded"]] = list(branded_lower,branded_lower,branded_upper,branded_upper)
  lcpars[["b_price"]] = list(price_lower,price_upper,price_lower,price_upper)
  
  
  
  classAlloc_settings = list(
    classes      = c(class =1,	class_2 =2,	class_3 =3,	class_4 =4), 
    avail        = 1
  )
  
  V = list()
  branded_c = (branded_a+branded_b)/2 + branded_c_delta
  price_c = (price_a+price_b)/2 + price_c_delta
  
  lower_branded=(branded_c-branded_a)/(branded_b-branded_a)
  upper_branded=(branded_b-branded_c)/(branded_b-branded_a)
  lower_price=(price_c-price_a)/(price_b-price_a)
  upper_price=(price_b-price_c)/(price_b-price_a)
  
  
  V[[1]] = log(lower_branded * lower_price)
  V[[2]] = log(lower_branded * upper_price)
  V[[3]] = log(upper_branded * lower_price)
  V[[4]] = log(upper_branded * upper_price)
  
  classAlloc_settings$utilities = V
  lcpars[["pi_values"]] = apollo_classAlloc(classAlloc_settings)
  
  return(lcpars)
}

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### Define settings for MNL model component that are generic across classes
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4),
    avail         = list(alt1=1, alt2=1, alt3=1, alt4=1),
    choiceVar     = new_choice
  )
  
  ### Loop over classes
  for(s in 1:4){
    ### Compute class-specific utilities
    V = list()
    V[["alt1"]] = ( branded[[s]]  * (brand_1=="Artemis"|brand_1=="Novum") 
                    + b_country_CH*(country_1=="Switzerland") + b_country_DK*(country_1=="Denmark") + b_country_USA*(country_1=="USA") 
                    + b_char_standard*(char_1=="standard") + b_char_fast*(char_1=="fast acting") + b_char_double*(char_1=="double strength") 
                    + b_risk*side_effects_1
                    + b_price[[s]]*price_1)
    V[["alt2"]] = ( branded[[s]]  * (brand_2=="Artemis"|brand_2=="Novum") 
                    + b_country_CH*(country_2=="Switzerland") + b_country_DK*(country_2=="Denmark") + b_country_USA*(country_2=="USA") 
                    + b_char_standard*(char_2=="standard") + b_char_fast*(char_2=="fast acting") + b_char_double*(char_2=="double strength") 
                    + b_risk*side_effects_2
                    + b_price[[s]]*price_2)
    V[["alt3"]] = ( b_unbranded *(brand_3=="BestValue"|brand_3=="Supermarket"|brand_3=="PainAway") 
                    + b_country_USA*(country_3=="USA") + b_country_IND*(country_3=="India") + b_country_RUS*(country_3=="Russia") + b_country_BRA*(country_3=="Brazil") 
                    + b_char_standard*(char_3=="standard") + b_char_fast*(char_3=="fast acting") 
                    + b_risk*side_effects_3
                    + b_price[[s]]*price_3 )
    V[["alt4"]] = ( b_unbranded*(brand_4=="BestValue"|brand_4=="Supermarket"|brand_4=="PainAway") 
                    + b_country_USA*(country_4=="USA") + b_country_IND*(country_4=="India") + b_country_RUS*(country_4=="Russia") + b_country_BRA*(country_4=="Brazil") 
                    + b_char_standard*(char_4=="standard") + b_char_fast*(char_4=="fast acting") 
                    + b_risk*side_effects_4
                    + b_price[[s]]*price_4 )
    
    mnl_settings$utilities = V

    ### Compute within-class choice probabilities using MNL model
    P[[paste0("Combination_",s)]] = apollo_mnl(mnl_settings, functionality)
    
    ### Take product across observation for same individual
    P[[paste0("Combination_",s)]] = apollo_panelProd(P[[paste0("Combination_",s)]], apollo_inputs ,functionality)
    
    ### Average across inter-individual draws within classes
    P[[paste0("Combination_",s)]] = apollo_avgInterDraws(P[[paste0("Combination_",s)]], apollo_inputs, functionality)
  }
  
  ### Compute latent class model probabilities
  lc_settings   = list(inClassProb = P, classProb=pi_values)
  P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #
### Optional starting values search
# apollo_beta=apollo_searchStart(apollo_beta, apollo_fixed,apollo_probabilities, apollo_inputs)
model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)


apollo_modelOutput(model)

#############
## LL
#############
apollo_beta = model[["estimate"]]
LL = log(apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"))
sum(LL) 

df1<-as.data.frame(1:2031)
df2<-as.data.frame(LL)
df<-cbind(df1,df2)
rm(df1,df2)
names(df)[1]<- "idd"
names(df)[2]<- "avgLL_AT"

setwd("J:/Model Averaging/Analysis/Empirical analysis/Tobacco choices")
write.csv(df,"SimDrugChoice2_AT_LL.csv", row.names=FALSE)
