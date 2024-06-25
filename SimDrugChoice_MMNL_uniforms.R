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

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "SimDrugChoice2_MMNL_uniforms",
  modelDescr      = "SimDrugChoice2_MMNL_uniforms",
  mixing    = TRUE, # TRUE for models that include random parameters
  nCores    = 3,    # Number of cores to use in estimation
  indivID         = "ID", 
  outputDirectory = "output"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
database = read.csv("Sim_Drug_Choice_2.csv")

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(
  lower_branded =3.250263,
  upper_branded =-4.054803,
  b_unbranded =0,
  lower_country_CH =-2.877799,
  upper_country_CH =0.562156,
  lower_country_DK =-0.430575,
  upper_country_DK =-1.989316,
  lower_country_USA =0,
  upper_country_USA =-2.413177,
  lower_country_IND =-3.830314,
  upper_country_IND =1.559533,
  lower_country_RUS =-1.945145,
  upper_country_RUS =-0.358573,
  lower_country_BRA =0.65113,
  upper_country_BRA =-3.101205,
  b_char_standard =0,
  lower_char_fast =0.394899,
  upper_char_fast =1.612473,
  lower_char_double =2.614585,
  upper_char_double =1.369035,
  lower_risk =-0.002732,
  upper_risk =-0.001042,
  lower_price =-1.210298,
  upper_price =-0.196472

  )

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_betlower_fixed = c() if none
apollo_fixed=names(apollo_beta)[apollo_beta==0|apollo_beta==1]


# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "mlhs", # use mlhs when using 5 or more distributions 
  interNDraws    = 500,
  interUnifDraws = c("draws_unif1_branded",
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
                     
                     "draws_unif2_branded",
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
  
  randcoeff[["b_branded"]]     = lower_branded     + (upper_branded     - lower_branded     ) * draws_unif1_branded 
  randcoeff[["b_country_CH"]]  = lower_country_CH  + (upper_country_CH  - lower_country_CH  ) * draws_unif1_CH  
  randcoeff[["b_country_DK"]]  = lower_country_DK  + (upper_country_DK  - lower_country_DK  ) * draws_unif1_DK      
  randcoeff[["b_country_USA"]] = lower_country_USA + (upper_country_USA - lower_country_USA ) * draws_unif1_USA    
  randcoeff[["b_country_IND"]] = lower_country_IND + (upper_country_IND - lower_country_IND ) * draws_unif1_IND    
  randcoeff[["b_country_RUS"]] = lower_country_RUS + (upper_country_RUS - lower_country_RUS ) * draws_unif1_RUS 
  randcoeff[["b_country_BRA"]] = lower_country_BRA + (upper_country_BRA - lower_country_BRA ) * draws_unif1_BRA
  randcoeff[["b_char_fast"]]   = lower_char_fast   + (upper_char_fast   - lower_char_fast   ) * draws_unif1_fast
  randcoeff[["b_char_double"]] = lower_char_double + (upper_char_double - lower_char_double ) * draws_unif1_double
  randcoeff[["b_risk"]]        = lower_risk        + (upper_risk        - lower_risk        ) * draws_unif1_risk
  randcoeff[["b_price"]]       = lower_price       + (upper_price       - lower_price       ) * draws_unif1_price
  
  return(randcoeff)
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
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  
  V[["alt1"]] = ( b_branded  * (brand_1=="Artemis"|brand_1=="Novum") 
                  + b_country_CH*(country_1=="Switzerland") + b_country_DK*(country_1=="Denmark") + b_country_USA*(country_1=="USA") 
                  + b_char_standard*(char_1=="standard") + b_char_fast*(char_1=="fast acting") + b_char_double*(char_1=="double strength") 
                  + b_risk*side_effects_1
                  + b_price*price_1)
  V[["alt2"]] = ( b_branded  * (brand_1=="Artemis"|brand_1=="Novum") 
                  + b_country_CH*(country_2=="Switzerland") + b_country_DK*(country_2=="Denmark") + b_country_USA*(country_2=="USA") 
                  + b_char_standard*(char_2=="standard") + b_char_fast*(char_2=="fast acting") + b_char_double*(char_2=="double strength") 
                  + b_risk*side_effects_2
                  + b_price*price_2)
  V[["alt3"]] = ( b_unbranded *(brand_3=="BestValue"|brand_3=="Supermarket"|brand_3=="PainAway") 
                  + b_country_USA*(country_3=="USA") + b_country_IND*(country_3=="India") + b_country_RUS*(country_3=="Russia") + b_country_BRA*(country_3=="Brazil") 
                  + b_char_standard*(char_3=="standard") + b_char_fast*(char_3=="fast acting") 
                  + b_risk*side_effects_3
                  + b_price*price_3 )
  V[["alt4"]] = ( b_unbranded*(brand_4=="BestValue"|brand_4=="Supermarket"|brand_4=="PainAway") 
                  + b_country_USA*(country_4=="USA") + b_country_IND*(country_4=="India") + b_country_RUS*(country_4=="Russia") + b_country_BRA*(country_4=="Brazil") 
                  + b_char_standard*(char_4=="standard") + b_char_fast*(char_4=="fast acting") 
                  + b_risk*side_effects_4
                  + b_price*price_4 )
  
  ### Compute probabilities for "best" choice using MNL model
  mnl_settings_best = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4),
    avail         = list(alt1=1, alt2=1, alt3=1, alt4=1),
    choiceVar     = new_choice,
    utilities     = V
  )
  P[["model"]] = apollo_mnl(mnl_settings_best, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #
### Estimate model
model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

s=list(printT1=TRUE)
apollo_modelOutput(model, modelOutput_settings=s)



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
names(df)[2]<- "avgLL_U"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
write.csv(df,"SimDrugChoice2_MMNL_uniforms_LL.csv", row.names=FALSE)


