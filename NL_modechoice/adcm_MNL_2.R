# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory
# rm(list = ls())

### Load Apollo library
library(apollo)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "MNL_2",
  modelDescr      = "MNL model with tradeoff parameter on mode choice RP data",
  indivID         = "NUMERIC", 
  outputDirectory = "output"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
### if data is to be loaded from a file (e.g. called data.csv), 
### the code would be: database = read.csv("data.csv",header=TRUE)
database = readxl::read_excel("Base Las Condes Centro.xls")

### Create a new variable of wage rate
database$w = database$ILM/(database$WS*4*60)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(asc_auto              = 0,
              asc_comp              = 0,
              asc_taxi              = 0,
              asc_metro             = 0,
              asc_bus               = 0,
              asc_autometro         = 0,
              asc_compmetro         = 0,
              asc_taximetro         = 0,
              asc_busmetro          = 0,

              #travel time
              b_tt_auto              = 0,
              b_tt_comp              = 0,
              b_tt_taxi              = 0,
              b_tt_metro             = 0,
              b_tt_bus               = 0,
              b_tt_autometro         = 0,
              b_tt_compmetro         = 0,
              b_tt_taximetro         = 0,
              b_tt_busmetro          = 0,

              #travel cost
              b_tc_auto              = 0,
              b_tc_comp              = 0,
              b_tc_taxi              = 0,
              b_tc_metro             = 0,
              b_tc_bus               = 0,
              b_tc_autometro         = 0,
              b_tc_compmetro         = 0,
              b_tc_taximetro         = 0,
              b_tc_busmetro          = 0,

              #access time
              b_acs_auto             = 0,
              b_acs_comp             = 0,
              b_acs_taxi             = 0,
              b_acs_metro            = 0,
              b_acs_bus              = 0,
              b_acs_autometro        = 0,
              b_acs_compmetro        = 0,
              b_acs_taximetro        = 0,
              b_acs_busmetro         = 0,

              #alighting time
              b_alt_taxi             = 0,
              b_alt_metro            = 0,
              b_alt_bus              = 0,
              
              alpha                  = 0)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_auto")

### Read in starting values for at least some parameters from existing model output file
# apollo_beta=apollo_readBeta(apollo_beta,apollo_fixed,"../output/MNL_SP_covariates",overwriteFixed=FALSE)

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
  
  ### Define alpha parameters
  a = exp(-alpha) /(1 + exp(-alpha))
  
  ### List of utilities: these must use the same names as in nl_settings, order is irrelevant
  V = list()
  V[["auto"     ]]  = asc_auto      + b_tt_auto      * TDV1 * (-w^a) + b_tc_auto      * CTOT1 * (-w^(a-1)) + b_acs_auto      * TCAM1
  V[["companion"]]  = asc_comp      + b_tt_comp      * TDV2 * (-w^a) + b_tc_comp      * CTOT2 * (-w^(a-1)) + b_acs_comp      * TCAM2
  V[["taxi"     ]]  = asc_taxi      + b_tt_taxi      * TDV3 * (-w^a) + b_tc_taxi      * CTOT3 * (-w^(a-1)) + b_acs_taxi      * TCAM3 + b_alt_taxi  * TESP3
  V[["metro"    ]]  = asc_metro     + b_tt_metro     * TDV4 * (-w^a) + b_tc_metro     * CTOT4 * (-w^(a-1)) + b_acs_metro     * TCAM4 + b_alt_metro * TESP4
  V[["bus"      ]]  = asc_bus       + b_tt_bus       * TDV5 * (-w^a) + b_tc_bus       * CTOT5 * (-w^(a-1)) + b_acs_bus       * TCAM5 + b_alt_bus   * TESP5
  V[["autometro"]]  = asc_autometro + b_tt_autometro * TDV6 * (-w^a) + b_tc_autometro * CTOT6 * (-w^(a-1)) + b_acs_autometro * TCAM6 + b_alt_metro * TESP4
  V[["compmetro"]]  = asc_compmetro + b_tt_compmetro * TDV7 * (-w^a) + b_tc_compmetro * CTOT7 * (-w^(a-1)) + b_acs_compmetro * TCAM7 + b_alt_metro * TESP4
  V[["taximetro"]]  = asc_taximetro + b_tt_taximetro * TDV8 * (-w^a) + b_tc_taximetro * CTOT8 * (-w^(a-1)) + b_acs_taximetro * TCAM8 + b_alt_metro * TESP4
  V[["busmetro" ]]  = asc_busmetro  + b_tt_busmetro  * TDV9 * (-w^a) + b_tc_busmetro  * CTOT9 * (-w^(a-1)) + b_acs_busmetro  * TCAM9 + b_alt_metro * TESP4
  
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(auto=1, companion=2, taxi=3, metro=4, bus=5, autometro=6, compmetro=7, taximetro=8, busmetro=9),
    avail        = list(auto=AVAIL1, companion=AVAIL2, taxi=AVAIL3, metro=AVAIL4, bus=AVAIL5, autometro=AVAIL6, compmetro=AVAIL7, taximetro=AVAIL8, busmetro=AVAIL9),
    choiceVar    = ICH,
    utilities     = V
  )
  
  ### Compute probabilities using MNL model
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  # P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model,modelOutput_settings = list(printT1=1))

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model,saveOutput_settings = list(printT1=1))

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

### Print outputs of additional diagnostics to new output file (remember to close file writing when complete)
apollo_sink()

# ----------------------------------------------------------------- #
#---- LR TEST AGAINST MNL MODEL                                  ----
# ----------------------------------------------------------------- #

apollo_lrTest("../output/MNL_1", model)

# ----------------------------------------------------------------- #
#---- switch off writing to file                                 ----
# ----------------------------------------------------------------- #

apollo_sink()
