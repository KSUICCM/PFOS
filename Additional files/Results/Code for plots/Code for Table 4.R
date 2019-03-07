################## loading all information and packages ######################
## Load libraries
library(mrgsolve)
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2) # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr) # Needed for the pipe %>% operator
library(FME) # Used to conduct local sensitivity analysis prior calibration and MCMC analysis
library(reshape) # melt function to reshape the table
library(minpack.lm) # For model calibration
library(truncnorm) # For truncated normal distribution in MCMC 
library(EnvStats) # For truncated normal distribution in MCMC
library(invgamma) # For MCMC simulation

## Input PBPK model (run the mrg code, and then saveRDS(RatPBPK.code, file = "ratPBPK.RDS"))
micePBPK.code     <- readRDS (file = "micePBPK.RDS")
monkeyPBPK.code   <- readRDS (file = "monkeyPBPK.RDS")
ratPBPK.code      <- readRDS (file = "ratPBPK.RDS")
humanPBPK.code    <- readRDS (file = "humanPBPK.RDS")

## Loading pbpk model (e.g., ratPBPK.code) to the object mod.rat (mcode is a function of loading model in mrgsolve)
mod.mouse         <- mcode ("micepbpk", micePBPK.code)
mod.rat           <- mcode ("ratpbpk", ratPBPK.code) 
mod.monkey        <- mcode ("monkeypbpk", monkeyPBPK.code)
mod.human         <- mcode ("humanpbpk", humanPBPK.code)


## Loading human, rat, mouse, monkey observed data
Human.obs         <- readRDS(file = "Human.obs.rds")
Rat.obs           <- readRDS(file = "Rat.obs.rds")
Mouse.obs         <- readRDS(file = "Mouse.obs.rds")
Monkey.obs        <- readRDS(file = "Monkey.obs.rds")

## Loading human, rat, mouse, monkey MCMC data (posterior parameter values from one of the four MCMC chains)
Human.MCMC        <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC          <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC        <- readRDS(file = "mouse.MCMC.rds")
Monkey.MCMC       <- readRDS(file = "Monkey.MCMC.rds")

## Population mean (u) of prior distribution (this part of the code (lines 39-47) is not needed in this file)
theta.Human       <- readRDS (file="theta.Human.Rds")
theta.Monkey      <- readRDS (file="theta.Monkey.Rds")
theta.Mouse       <- readRDS (file="theta.Mouse.Rds")
theta.Rat         <- readRDS (file="theta.Rat.Rds")

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)


######################### Table 4: Comparision of the POD-derived HED with EPA results ####################
## Predition function

## Prediction function for monkeys, rats, and humans
pred.AUC <- function(pars.monkey,pars.rat,pars.human) {
  
  pars.monkey %<>% lapply(exp)
  names(pars.monkey) <- names(pars.monkey)
  pars.monkey        <- pars.monkey [-which_sig ] # pars.monkey, pars.rat, pars.human are defined in lines 182-184.
  
  ## Define the three exposure scenarios:
  ## Monkeys were administered 0.15 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  BW.monkey.A           = 3.9                                  ## Monkey body weight

  tinterval.monkey      = 24                                   ## Time interval
  TDoses.monkey.rep     = 182                                 ## Total number of doses, daily dosing for 182 days
  
  PDOSEoral.monkey.A    = 0.15                                 ## mg/kg; BW Oral dose
  DOSEoral.monkey.rep.A = PDOSEoral.monkey.A*BW.monkey.A       ## mg; amount of oral dose

  ex.monkey.rep.A <- ev(ID = 1, amt= DOSEoral.monkey.rep.A, ii = tinterval.monkey, addl = TDoses.monkey.rep-1, cmt="AST", replicate = FALSE)

  ## Set up the exposure time
  tsamp.monkey.rep = tgrid(0,tinterval.monkey*(TDoses.monkey.rep-1)+24*1,24)  ## 182 days and simulated for 24*565 hours after dosing, output interval is 24 hours
  
  out.monkey <- 
    mod.monkey  %>% 
    param(pars.monkey) %>%
    Req (Plasma, AUC_CA,AUC_CL)%>% # Plasma and AUC_CA are defined in the mrg model code.
    update(atol = 1E-10, maxsteps = 10000) %>% # simulation time step, max steps per hour
    mrgsim_d (data = ex.monkey.rep.A, tgrid = tsamp.monkey.rep) 
  
  outdf.monkey  <- cbind.data.frame( # combine output data into data frame for plotting in ggplot
                   Time   = out.monkey$time/24,                                # Sample the data every 24 hours
                   CA     = out.monkey$Plasma,                                 ## ug/ml, PFOS levels in plasma
                   AUC.CA = out.monkey$AUC_CA,
                   AUC.CL = out.monkey$AUC_CL)
  
  ## Define the exposure scenario for rats: oral exposure to 0.34 mg/kg 
  
  pars.rat %<>% lapply(exp)
  names(pars.rat) <- names(pars.rat)
  pars.rat        <- pars.rat [-which_sig ]
  

  BW.rat                 = 0.3                                       # Rat body weight
  tinterval.rat          = 24                                        # Time interval
  TDoses.rat             = 98                                        # Dose times (total number of dosing)
  PDOSEoral.rat          = 0.34                                      # This value is from Table 4-8 (Page 4-14) in the 2016 EPA Report, which was based on Seacat et al. 2003
  DOSEoral.rat           = PDOSEoral.rat*BW.rat                      # amount of oral dose
  ex.rat                 <- ev (ID = 1, amt= DOSEoral.rat, 
                                ii = tinterval.rat, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  tsamp.rat.rep          = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*1,24)   ## Siumuation 24*98 + 24*100 hours (198 days)
  
  ## Prediction of AUC
  out.rat <- 
    mod.rat %>% 
    param(pars.rat) %>%
    Req (Plasma, AUC_CA, AUC_CL)%>%
    update(atol = 1E-10, maxsteps = 10000) %>%
    mrgsim_d (data = ex.rat, tgrid = tsamp.rat.rep)
  
  outdf.rat  <- cbind.data.frame(
                Time   = out.rat$time/24,
                CA     = out.rat$Plasma,                                 ## ug/ml,PFOS levels in plasma
                AUC.CA = out.rat$AUC_CA,
                AUC.CL = out.rat$AUC_CL)
  
  
  ## Human scenario A: oral exposue to 0.15 mg/kg-d
  pars.human %<>% lapply(exp)
  names(pars.human) <- names(pars.human)
  pars.human        <- pars.human [-which_sig ]
  
  ###
  BW.human             = 82.3                                 ## kg

  ## 
  tinterval.human      = 24                                   ## Time interval
  TDoses.human.A       = 365*25                               ## Dose time for 25 years to reach the steady state
  TDoses.human.B       = 365*25                               ## Dose time for 25 years to reach the steady state
  
  PDOSEoral.human.A    = 0.15                                  ## mg/kg; BW Oral dose
  PDOSEoral.human.B    = 0.34                                  ## mg/kg; BW Oral dose
  
  DOSEoral.human.A = PDOSEoral.human.A* BW.human              ## mg; amount of oral dose
  DOSEoral.human.B = PDOSEoral.human.B* BW.human              ## mg; amount of oral dose
  
  ex.human.A <- ev(ID=1, amt= DOSEoral.human.A, ii = tinterval.human, addl = TDoses.human.A-1, cmt="AST", replicate = FALSE)
  ex.human.B <- ev(ID=1, amt= DOSEoral.human.B, ii = tinterval.human, addl = TDoses.human.B-1, cmt="AST", replicate = FALSE)
  
  ## set up the exposure time
  tsamp.human.A = tgrid(0,tinterval.human*(TDoses.human.A-1)+24*1,24)  ## 25 years and simulated for 24*365 hours (1 year) after dosing
  tsamp.human.B = tgrid(0,tinterval.human*(TDoses.human.B-1)+24*1,24)  ## 25 years and simulated for 24*365 hours (1 year) after dosing
  
  
  out.human.A <- 
    mod.human %>% 
    param(pars.human) %>%
    Req (Plasma, AUC_CA, AUC_CL)%>%
    update(atol = 1E-10,rtol= 1e-10, maxsteps = 50000) %>%
    mrgsim_d(data = ex.human.A, tgrid = tsamp.human.A)
  
  outdf.human.A <- cbind.data.frame(
                Time   = out.human.A$time/24,
                CA     = out.human.A$Plasma,                                 ## ug/ml,PFOS levels in plasma
                AUC.CA = out.human.A$AUC_CA,
                AUC.CL = out.human.A$AUC_CL)
  
  out.human.B <- 
    mod.human %>% 
    param(pars.human) %>%
    Req (Plasma, AUC_CA, AUC_CL)%>%
    update(atol = 1E-10,rtol= 1e-10, maxsteps = 50000) %>%
    mrgsim_d(data = ex.human.B, tgrid = tsamp.human.B)
  
  outdf.human.B <- cbind.data.frame(
    Time   = out.human.B$time/24,
    CA     = out.human.B$Plasma,                                 ## ug/ml,PFOS levels in plasma
    AUC.CA = out.human.B$AUC_CA,
    AUC.CL = out.human.B$AUC_CL)
  
  
  return(list("monkey"  = outdf.monkey,
              "rat"     = outdf.rat, # save the outdf.rat result to the object of rat
              "human.m" = outdf.human.A,
              "human.r" = outdf.human.B))
  
}

## Probablity of AUC
AUC = matrix(nrow = 5000, ncol = 8) # Create an empty matrix to store the simulation data.

for (i in 1:5000){
  
  j = i*10 
  pars.monkey = Monkey.MCMC$pars [j,]
  pars.rat    = Rat.MCMC$pars [j,]
  pars.human  = Human.MCMC$pars [j,]
  pars.human [2] = log(exp(pars.human [2])/1000) # km of human is ug/L and transfered to unit of mg/L corresponding to the unit of animals
  pars.human [4] = log(exp(pars.human [4])/1000)
  
  MC.HED      <- pred.AUC(pars.monkey,pars.rat,pars.human) # Generate the results of 5000 Plasma and AUC_CA. Plasma is included here in case we want to use Cmax for deriving POD
  AUC   [i,1] <- MC.HED$monkey [MC.HED$monkey$Time == 182, ]$AUC.CA     ## AUC of plasma in monkeys based on NOAEL. To obtain AUC.CA for each iteration for a total of 5000 iterations
  AUC   [i,2] <- MC.HED$rat    [MC.HED$rat$Time == 98, ]$AUC.CA         ## AUC of plasma in rats based on NOAEL 
  AUC   [i,3] <- MC.HED$human.m[MC.HED$human.m$Time == 365*25, ]$AUC.CA    ## AUC of plasma in humans based on monkey NOAEL 
  AUC   [i,4] <- MC.HED$human.r[MC.HED$human.r$Time == 365*25, ]$AUC.CA     ## AUC of palsma in humans based on rat NOAEL
  AUC   [i,5] <- MC.HED$monkey [MC.HED$monkey$Time == 182, ]$AUC.CL     ## AUC of liver in monkeys based on NOAEL. To obtain AUC.CL for each iteration for a total of 5000 iterations
  AUC   [i,6] <- MC.HED$rat    [MC.HED$rat$Time == 98, ]$AUC.CL         ## AUC of liver in rats based on NOAEL 
  AUC   [i,7] <- MC.HED$human.m[MC.HED$human.m$Time == 365*25, ]$AUC.CL    ## AUC of liver in humans based on monkey NOAEL 
  AUC   [i,8] <- MC.HED$human.r[MC.HED$human.r$Time == 365*25, ]$AUC.CL     ## AUC of liver in humans based on rat NOAEL
  cat("iteration = ", i , "\n") # caption the result of each iteration per line
}


######################################################
colnames(AUC) = c("AUC.CA.Monkey","AUC.CA.rat","AUC.CA.human.m","AUC.CA.human.r",
                  "AUC.CL.Monkey","AUC.CL.rat","AUC.CL.human.m","AUC.CL.human.r")

AUC  = as.data.frame(AUC)  

## Estimated the average serum concentraiton (ASC) based on the rat and moneky NOAEL
AUC$Avg.CA.monkey    = AUC$AUC.CA.Monkey/(182*24)      # Monkey ASC: AUC normalized by expsoure duration (182 days or 182*24 hours) 
AUC$Avg.CA.rat       = AUC$AUC.CA.rat/(98*24)          # Rat ASC: AUC normalized by expsoure duration (98 days or 98*24 hours) 
AUC$Avg.CA.human.m   = AUC$AUC.CA.human.m/(365*25*24)  # Human ASC: AUC normalized by the duration reaching the steady state (Based on the simulation results, it reached steady state at ~25 years.)   
AUC$Avg.CA.human.r   = AUC$AUC.CA.human.r/(365*25*24)

## Estimated the average liver concentration (ALC)
AUC$Avg.CL.monkey    = AUC$AUC.CL.Monkey/(182*24)     
AUC$Avg.CL.rat       = AUC$AUC.CL.rat/(98*24)          
AUC$Avg.CL.human.m   = AUC$AUC.CL.human.m/(365*25*24) 
AUC$Avg.CL.human.r   = AUC$AUC.CL.human.r/(365*25*24)


## Calculated the HED 
## This study method
AUC$HED.CA.human.m      = ((AUC$Avg.CA.monkey)/(AUC$Avg.CA.human.m))*0.15 # monkey NOAEL: 0.15 (Seacat et al., 2002) POD_human = NOAEL_human = NOAEL_animal * (ASC_animal/ASC_human)
AUC$HED.CA.human.r      = ((AUC$Avg.CA.rat)/(AUC$Avg.CA.human.r))*0.34    # rat NoAEL: 0.34 (Seacat et al., 2003)
AUC$HED.CL.human.m      = ((AUC$Avg.CL.monkey)/(AUC$Avg.CL.human.m))*0.15 # monkey NOAEL: 0.15 (Seacat et al., 2002) POD_human = NOAEL_human = NOAEL_animal * (ASC_animal/ASC_human)
AUC$HED.CL.human.r      = ((AUC$Avg.CL.rat)/(AUC$Avg.CL.human.r))*0.34    # rat NoAEL: 0.34 (Seacat et al., 2003)


## Estimated the median (95% CI) of ASC
# Monkey
ASC.m.range <- cbind.data.frame (
  Median.CA  = quantile(AUC$Avg.CA.monkey  , probs = 0.5, names = FALSE,na.rm=T),
  lower.CA   = quantile(AUC$Avg.CA.monkey  , probs = 0.025, names = FALSE,na.rm=T),
  upper.CA   = quantile(AUC$Avg.CA.monkey  , probs = 0.975, names = FALSE,na.rm=T),
  Median.CL  = quantile(AUC$Avg.CL.monkey  , probs = 0.5, names = FALSE,na.rm=T),
  lower.CL   = quantile(AUC$Avg.CL.monkey  , probs = 0.025, names = FALSE,na.rm=T),
  upper.CL   = quantile(AUC$Avg.CL.monkey  , probs = 0.975, names = FALSE,na.rm=T)
)

# Rat
ASC.r.range <- cbind.data.frame (
  Median.CA    = quantile(AUC$Avg.CA.rat, probs = 0.5, names = FALSE,na.rm=T),
  lower.CA     = quantile(AUC$Avg.CA.rat, probs = 0.025, names = FALSE,na.rm=T),
  upper.CA     = quantile(AUC$Avg.CA.rat, probs = 0.975, names = FALSE,na.rm=T),
  Median.CL    = quantile(AUC$Avg.CL.rat, probs = 0.5, names = FALSE,na.rm=T),
  lower.CL     = quantile(AUC$Avg.CL.rat, probs = 0.025, names = FALSE,na.rm=T),
  upper.CL     = quantile(AUC$Avg.CL.rat, probs = 0.975, names = FALSE,na.rm=T)
)

## Estimated HED based on the AUC Method
# Monkey
HED.human.m.range <- cbind.data.frame (
Median.CA  = quantile(AUC$HED.CA.human.m  , probs = 0.5, names = FALSE,na.rm=T),
lower.CA   = quantile(AUC$HED.CA.human.m  , probs = 0.025, names = FALSE,na.rm=T),
upper.CA   = quantile(AUC$HED.CA.human.m  , probs = 0.975, names = FALSE,na.rm=T),
Median.CL  = quantile(AUC$HED.CL.human.m  , probs = 0.5, names = FALSE,na.rm=T),
lower.CL   = quantile(AUC$HED.CL.human.m  , probs = 0.025, names = FALSE,na.rm=T),
upper.CL   = quantile(AUC$HED.CL.human.m  , probs = 0.975, names = FALSE,na.rm=T)
)

# Rat
HED.human.r.range <- cbind.data.frame (
Median.CA  = quantile(AUC$HED.CA.human.r , probs = 0.5, names = FALSE,na.rm=T),
lower.CA   = quantile(AUC$HED.CA.human.r , probs = 0.025, names = FALSE,na.rm=T),
upper.CA   = quantile(AUC$HED.CA.human.r , probs = 0.975, names = FALSE,na.rm=T),
Median.CL  = quantile(AUC$HED.CL.human.r , probs = 0.5, names = FALSE,na.rm=T),
lower.CL   = quantile(AUC$HED.CL.human.r , probs = 0.025, names = FALSE,na.rm=T),
upper.CL   = quantile(AUC$HED.CL.human.r , probs = 0.975, names = FALSE,na.rm=T)
)

