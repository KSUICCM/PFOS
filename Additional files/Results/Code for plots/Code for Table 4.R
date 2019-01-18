################## loading all information and packages ######################
## Load libraries
library(mrgsolve)
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2) # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr) # Needed for the pipe %>% operator
library(FME) # Used to conduct local sensitivity analysis prior calibration and MCMC analysis
library(reshape) # melt function to reshape the table
library(minpack.lm)
library(truncnorm)   
library(EnvStats)
library(invgamma)

## input PBPK model (run the mrg code, and then saveRDS(RatPBPK.code, file = "ratPBPK.RDS"))
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


######################### Table 4: Comparision the HED POD with EPA ####################
## Predition function

## Prediction function for monkey, rat, and human
pred.AUC <- function(pars.monkey,pars.rat,pars.human) {
  
  pars.monkey %<>% lapply(exp)
  names(pars.monkey) <- names(pars.monkey)
  pars.monkey        <- pars.monkey [-which_sig ]
  
  ## Define the three exposure scenario:
  ## Monkeys were administered 0.15 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  BW.monkey.A           = 3.9                                  ## Monkey body weight

  tinterval.monkey      = 24                                   ## Time interval
  TDoses.monkey.rep     = 182                                  ## Total number of doses, daily dosing for 182 days
  
  PDOSEoral.monkey.A    = 0.15                                 ## mg/kg; BW Oral dose
  DOSEoral.monkey.rep.A = PDOSEoral.monkey.A*BW.monkey.A       ## mg; amount of oral dose

  ex.monkey.rep.A <- ev(ID = 1, amt= DOSEoral.monkey.rep.A, ii = tinterval.monkey, addl = TDoses.monkey.rep-1, cmt="AST", replicate = FALSE)

  ## set up the exposure time
  tsamp.monkey.rep = tgrid(0,tinterval.monkey*(TDoses.monkey.rep-1)+24*565,24)  ## 182 days and simulated for 24*565 hours after dosing, output interval is 24 hours
  
  out.monkey <- 
    mod.monkey  %>% 
    param(pars.monkey) %>%
    Req (Plasma, AUC_CA)%>%
    update(atol = 1E-10, maxsteps = 10000) %>% # simulation time step, max steps per hour
    mrgsim_d (data = ex.monkey.rep.A, tgrid = tsamp.monkey.rep) 
  
  outdf.monkey  <- cbind.data.frame( # combine output data into data frame for plotting in ggplot
                   Time   = out.monkey$time/24,
                   CA     = out.monkey$Plasma,                                 ## ug/ml, PFOS levels in plasma
                   AUC.CA = out.monkey$AUC_CA)
  
  ## Define the exposure scenario for rats: oral exposure to 0.34  
  
  pars.rat %<>% lapply(exp)
  names(pars.rat) <- names(pars.rat)
  pars.rat        <- pars.rat [-which_sig ]
  

  BW.rat                 = 0.3                                       # Rat body weight
  tinterval.rat          = 24                                        # Time interval
  TDoses.rat             = 98                                         # Dose times (total number of dosing)
  PDOSEoral.rat          = 0.34                                       # Seacat et al. 2003
  DOSEoral.rat           = PDOSEoral.rat*BW.rat                        # amount of oral dose
  ex.rat                 <- ev (ID = 1, amt= DOSEoral.rat, 
                                ii = tinterval.rat, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  tsamp.rat.rep          = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*100,1)   ## Siumuation 24*98 + 24*100 hours (198 days)
  
  ## Prediction AUC
  out.rat <- 
    mod.rat %>% 
    param(pars.rat) %>%
    Req (Plasma, AUC_CA)%>%
    update(atol = 1E-10, maxsteps = 10000) %>%
    mrgsim_d (data = ex.rat, tgrid = tsamp.rat.rep)
  
  outdf.rat  <- cbind.data.frame(
                Time   = out.rat$time/24,
                CA     = out.rat$Plasma,                                 ## ug/ml,PFOS levels in plasma
                AUC.CA = out.rat$AUC_CA)
  
  
  ## Human scenario A: oral exposue to 0.15 mg/kg-d
  pars.human %<>% lapply(exp)
  names(pars.human) <- names(pars.human)
  pars.human        <- pars.human [-which_sig ]
  
  ###
  BW.human             = 82.3

  ## 
  tinterval.human      = 24                                   ## Time interval
  TDoses.human.A       = 182                                  ## Dose times 182 days
  TDoses.human.B       = 98                                   ## Dose times 98 days
  
  PDOSEoral.human.A    = 150                                  ## ug/kg; BW Oral dose
  PDOSEoral.human.B    = 340                                  ## ug/kg; BW Oral dose
  
  DOSEoral.human.A = PDOSEoral.human.A* BW.human              ## ug; amount of oral dose
  DOSEoral.human.B = PDOSEoral.human.B* BW.human              ## ug; amount of oral dose
  
  ex.human.A <- ev(ID=1, amt= DOSEoral.human.A, ii = tinterval.human, addl = TDoses.human.A-1, cmt="AST", replicate = FALSE)
  ex.human.B <- ev(ID=1, amt= DOSEoral.human.B, ii = tinterval.human, addl = TDoses.human.B-1, cmt="AST", replicate = FALSE)
  
  ## set up the exposure time
  tsamp.human.A = tgrid(0,tinterval.human*(TDoses.human.A-1)+24*100,24)  ## 182 days and simulated for 24*100 hours after dosing
  tsamp.human.B = tgrid(0,tinterval.human*(TDoses.human.B-1)+24*100,24)  ## 98 days and simulated for 24*100 hours after dosing
  
  
  out.human.A <- 
    mod.human %>% 
    param(pars.human) %>%
    Req (Plasma, AUC_CA)%>%
    update(atol = 1E-10,rtol= 1e-10, maxsteps = 50000) %>%
    mrgsim_d(data = ex.human.A, tgrid = tsamp.human.A)
  
  outdf.human.A <- cbind.data.frame(
                Time   = out.human.A$time/24,
                CA     = out.human.A$Plasma,                                 ## ug/ml,PFOS levels in plasma
                AUC.CA = out.human.A$AUC_CA)
  
  out.human.B <- 
    mod.human %>% 
    param(pars.human) %>%
    Req (Plasma, AUC_CA)%>%
    update(atol = 1E-10,rtol= 1e-10, maxsteps = 50000) %>%
    mrgsim_d(data = ex.human.B, tgrid = tsamp.human.B)
  
  outdf.human.B <- cbind.data.frame(
    Time   = out.human.B$time/24,
    CA     = out.human.B$Plasma,                                 ## ug/ml,PFOS levels in plasma
    AUC.CA = out.human.B$AUC_CA)
  
  
  
  
  return(list("monkey"  = outdf.monkey,
              "rat"     = outdf.rat, # save the outdf.rat result to the object of rat
              "human.m" = outdf.human.A,
              "human.r" = outdf.human.B))
  
}

## Probablity of AUC
AUC = matrix(nrow = 5000, ncol = 4)

for (i in 1:5000){
  
  j = i*10 
  pars.monkey = Monkey.MCMC$pars [j,]
  pars.rat    = Rat.MCMC$pars [j,]
  pars.human  = Human.MCMC$pars [j,]
  
  MC.AUC     <- pred.AUC(pars.monkey,pars.rat,pars.human)
  AUC  [i,1] <- MC.AUC$monkey [MC.AUC$monkey$Time == 182, ]$AUC.CA     ## AUC of plasma in monkey based on NOAEL 
  AUC  [i,2] <- MC.AUC$rat    [MC.AUC$rat$Time == 98, ]$AUC.CA         ## AUC of plasma in RAT based on NOAEL 
  AUC  [i,3] <- MC.AUC$human.m[MC.AUC$human.m$Time == 182, ]$AUC.CA    ## AUC of plasma in Human based on monkey NOAEL 
  AUC  [i,4] <- MC.AUC$human.r[MC.AUC$human.r$Time == 98, ]$AUC.CA     ## AUC of palsma in Human based on rat NOAEL
  cat("iteration = ", i , "\n") # caption the result of each iteration per line
}


####################
colnames(AUC) = c("AUC.Monkey","AUC.rat",
                  "AUC.human.m","AUC.human.r")

AUC = as.data.frame(AUC)  

## Estimated the Average serum concentraiton based on rat and moneky NOALE
AUC$Avg.CA.monkey = AUC$AUC.Monkey/(182*24)  # will create a new column of Avg.CA.monkey in AUC data frame
AUC$Avg.CA.rat    = AUC$AUC.rat/(98*24)
AUC$Avg.human.m   = AUC$AUC.human.m/(182*24*1000)
AUC$Avg.human.r   = AUC$AUC.human.r/(98*24*1000)

## Calculated the HED 
AUC$HED.monkey      = AUC$Avg.CA.monkey*0.000081 
AUC$HED.rat         = AUC$Avg.CA.rat*0.000081 
#AUC$HED.human.m     = (AUC$Avg.human.m)/(AUC$Avg.CA.monkey)*0.15
#AUC$HED.human.r     = (AUC$Avg.human.r)/(AUC$Avg.CA.rat)*0.34 



## Estimated the median (95% CI)
Avg.human.m.range <- cbind.data.frame (
  Median  = quantile(AUC$Avg.CA.monkey  , probs = 0.5, names = FALSE,na.rm=T),
  lower   = quantile(AUC$Avg.CA.monkey  , probs = 0.025, names = FALSE,na.rm=T),
  upper   = quantile(AUC$Avg.CA.monkey  , probs = 0.975, names = FALSE,na.rm=T)
)

Avg.human.r.range <- cbind.data.frame (
  Median  = quantile(AUC$Avg.CA.rat   , probs = 0.5, names = FALSE,na.rm=T),
  lower   = quantile(AUC$Avg.CA.rat   , probs = 0.025, names = FALSE,na.rm=T),
  upper   = quantile(AUC$Avg.CA.rat   , probs = 0.975, names = FALSE,na.rm=T)
)

# HED.human.m.range <- cbind.data.frame (
# Median  = quantile(AUC$HED.human.m  , probs = 0.5, names = FALSE,na.rm=T),
# lower   = quantile(AUC$HED.human.m  , probs = 0.025, names = FALSE,na.rm=T),
# upper   = quantile(AUC$HED.human.m  , probs = 0.975, names = FALSE,na.rm=T)
# )
# 
# HED.human.r.range <- cbind.data.frame (
# Median  = quantile(AUC$HED.human.r , probs = 0.5, names = FALSE,na.rm=T),
# lower   = quantile(AUC$HED.human.r , probs = 0.025, names = FALSE,na.rm=T),
# upper   = quantile(AUC$HED.human.r , probs = 0.975, names = FALSE,na.rm=T)
# )

HED.monkey.range <- cbind.data.frame (
  Median  = quantile(AUC$HED.monkey   , probs = 0.5, names = FALSE,na.rm=T),
  lower   = quantile(AUC$HED.monkey   , probs = 0.025, names = FALSE,na.rm=T),
  upper   = quantile(AUC$HED.monkey   , probs = 0.975, names = FALSE,na.rm=T)
)

HED.rat.range <- cbind.data.frame (
  Median  = quantile(AUC$HED.rat , probs = 0.5, names = FALSE,na.rm=T),
  lower   = quantile(AUC$HED.rat , probs = 0.025, names = FALSE,na.rm=T),
  upper   = quantile(AUC$HED.rat , probs = 0.975, names = FALSE,na.rm=T)
)


