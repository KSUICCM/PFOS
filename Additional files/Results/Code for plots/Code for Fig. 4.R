########################################## Fig. 4 ################################# 
# Fig. 4a Global fitting analysis (to generate the goodness of fit plot)          #
# Fig. 4b Predicted-to-observed ratio residuals versus model prediction plot      #
###################################################################################
## loading R packages
## Load libraries
library(mrgsolve) # Needed to run the main PBPK code
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)  # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr)    # Needed for the pipe %>% operator
library(reshape)  # melt function to reshape the table, reshape2 is version 2 of reshape. reshape is more stable when using melt function. We decided to use reshape for all.

## Input PBPK model
micePBPK.code     <- readRDS (file = "micePBPK.RDS")
monkeyPBPK.code   <- readRDS (file = "monkeyPBPK.RDS")
ratPBPK.code      <- readRDS (file = "ratPBPK.RDS")
humanPBPK.code    <- readRDS (file = "humanPBPK.RDS")

## Loading pbpk model 
mod.mouse         <- mcode ("micepbpk", micePBPK.code)
mod.rat           <- mcode ("ratpbpk", ratPBPK.code)
mod.monkey        <- mcode ("monkeypbpk", monkeyPBPK.code)
mod.human         <- mcode ("humanpbpk", humanPBPK.code)


## Loading human, rat, mouse, monkey observed data
Human.obs         <- readRDS(file = "Human.obs.rds")
Rat.obs           <- readRDS(file = "Rat.obs.rds")
Mouse.obs         <- readRDS(file = "Mouse.obs.rds")
Monkey.obs        <- readRDS(file = "Monkey.obs.rds")

## Loading human, rat, mouse, monkey MCMC data
Human.MCMC        <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC          <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC        <- readRDS(file = "Mouse.MCMC.rds")
Monkey.MCMC       <- readRDS(file = "Monkey.MCMC.rds")

## Population mean (u)
theta.Human  = readRDS (file="theta.Human.Rds")
theta.Monkey = readRDS (file="theta.Monkey.Rds")
theta.Mouse  = readRDS (file="theta.Mouse.Rds")
theta.Rat    = readRDS (file="theta.Rat.Rds")

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)

############################# prediction function ############################
############################################# Mouse observed data #####################################################                                                                                                                 #      
# A1: CD1 mice oral single dose of 1  mg/kg; plasma       Chang et al. (2012)                                         #
# A2: CD1 mice oral single dose of 1  mg/kg; liver                                                                    #
# A3: CD1 mice oral single dose of 1  mg/kg; kidney                                                                   #
# A4: CD1 mice oral single dose of 1  mg/kg; urine                                                                    # 
# B1: CD1 mice oral single dose of 20 mg/kg; plasma                                                                   #
# B2: CD1 mice oral single dose of 20 mg/kg; liver                                                                    #
# B3: CD1 mice oral single dose of 20 mg/kg; kidney                                                                   #
# B4: CD1 mice oral single dose of 20 mg/kg; urine                                                                    #
#######################################################################################################################

## Prediction function for mouse
pred.mouse <- function(pars, pred=FALSE) { # if pred=true, then generate predicted data for the entire time course; if false, then only generate predicted data at observed time points.
  
  ## Define the exposure scenario for oral exposue to 1 and 20 mg/kg-d
  BW.mouse             = 0.025                                   ## Mouse body weight
  tinterval.mouse      = 24                                      ## Time interval
  TDoses.mouse         = 1                                       ## Total number of dosing
  PDOSEoral.mouse.A    = 1                                       ## 1  mg/kg-d, Single oral dose
  PDOSEoral.mouse.B    = 20                                      ## 20 mg/kg-d, Single oral dose
  
  DOSEoral.mouse.A     = PDOSEoral.mouse.A*BW.mouse              ## Amount of oral dose
  DOSEoral.mouse.B     = PDOSEoral.mouse.B*BW.mouse              ## Amount of oral dose
  
  ex.mouse.A           <- ev (ID = 1, amt= DOSEoral.mouse.A, 
                              ii = tinterval.mouse, addl=TDoses.mouse-1, cmt="AST", replicate = FALSE)
  ex.mouse.B           <- ev (ID = 1, amt= DOSEoral.mouse.B, 
                              ii = tinterval.mouse, addl=TDoses.mouse-1, cmt="AST", replicate = FALSE)
  
  
  ## set up the exposure time
  tsamp.mouse          = tgrid(0,tinterval.mouse*(TDoses.mouse-1)+24*200,1)     ## Simulation for 24*200 hours (200 days), sample the data every hour
  
  
  ## Get out of log domain
  pars %<>% lapply(exp) # apply exp function to the object pars, and then assign the new values to the object pars
  names(pars) <- names(pars)
  
  ## Get a prediction
  out.A <- 
    mod.mouse %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>% # Plasma, Liver, and Kidney are defined in the mrg code (lines 182-184 for the rat)
    update(atol = 1E-8,maxsteps = 5000) %>% # simulation time step, max steps per hour
    mrgsim_d (data = ex.mouse.A, tgrid = tsamp.mouse) 
  
  outdf.A <- cbind.data.frame(Time   = out.A$time/24, 
                              CA     = out.A$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.A$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.A$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine = (out.A$Aurine/DOSEoral.mouse.A)*100)     ## %Dose, PFOS in urine
  out.B <- 
    mod.mouse %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>% # Plasma, Liver, and Kidney are defined in the mrg code (lines 182-184 for the rat)
    update(atol = 1E-8, maxsteps = 5000) %>% # simulation time step, max steps per hour
    mrgsim_d (data = ex.mouse.B, tgrid = tsamp.mouse)
  
  outdf.B <- cbind.data.frame(Time   = out.B$time/24, 
                              CA     = out.B$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.B$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.B$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine = (out.B$Aurine/DOSEoral.mouse.B)*100)     ## %Dose, PFOS in urine
  
  if (pred) return (list("outdf.A"= outdf.A, "outdf.B" = outdf.B)) # if shows the "~.." error code, just delete it; if pred=true, then run this line, generate entire time course predictions
  # if false, then run the lines below
  
  
  
  outdf.A = outdf.A[which(outdf.A$Time %in% Mouse.obs$A1$Time),] # "%in%" overlapped the observed and predicted data at the observed time points into the same object
  outdf.B = outdf.B[which(outdf.B$Time %in% Mouse.obs$B1$Time),]
  
  return (list("outdf.A"= outdf.A, "outdf.B" = outdf.B))
  
}



############################################# Rat data ################################################################
# Rat.obs:                                                                                                            #    
# A1: SD rat oral single dose to 4.2 mg/kg, matrix: plasma, data from Luccisano et al. (2012); Chang et al. (2012)    #                  
# A2: SD rat oral single dose to 2   mg/kg, matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #
# A3: SD rat oral single dose to 2   mg/kg, matrix: plasma, data from kim et al. (2016)                               #
# A4: SD rat iv exposure to 4.2 mg/kg,      matrix: urine, data from Johnson et al., 1979;                            #
# B1: SD rat oral single dose to 15 mg/kg,  matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #
# B2: SD rat oral single dose to 15 mg/kg,  matrix: urine, data from Chang et al. (2012)                              #
# B3: SD rat iv single dose to 2 mg/kg,     matrix: plasma, data from kim et al. (2016)                               #
# B4: SD rat oral daily dose to 1 mg/kg for 28 days,matrix: plasma,data from 3M unpublished data                      #
#######################################################################################################################

## Preiction function
pred.rat <- function(pars,pred=FALSE) {
  
  ## Exposure scenario for oral exposue to 4.2 mg/kg
  
  BW.rat                 = 0.3                                       # Rat body weight
  tinterval.rat          = 24                                        # Time interval
  TDoses.rat             = 1                                         # Dose times
  PDOSEoral.rat.A        = 4.2                                       # Single oral dose from Luccisano et al. (2012)
  DOSEoral.rat.A         = PDOSEoral.rat.A*BW.rat                    # Amount of oral dose
  ex.rat.A               <- ev (ID = 1, amt= DOSEoral.rat.A, 
                                ii = tinterval.rat, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  PDOSEoral.rat.B        = 2                                         ## Single oral dose from Luccisano et al. (2012)
  DOSEoral.rat.B         = PDOSEoral.rat.B*BW.rat                    ## Amount of oral dose
  ex.rat.B               <- ev (ID=1, amt= DOSEoral.rat.B, 
                                ii=tinterval.rat, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  
  ## Exposure scenario for oral exposue to 15 mg/kg
  PDOSEoral.rat.C        = 15                                        ## mg/kg, single oral dose                           
  DOSEoral.rat.C         = PDOSEoral.rat.C*BW.rat                    ## mg, amount of oral dose                   
  ex.rat.C               <- ev(ID=1,amt= DOSEoral.rat.C, 
                               ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  
  ## Exposure scenario for oral daily exposue to 1 mg/kg-d for 28 days
  tinterval.rep          = 24                                                    ## hr,    time interval
  TDoses.rep             = 28                                                    ## times, Dose times
  
  PDOSEoral.rat.D        = 1                                                     ## mg/kg, Oral dose
  DOSEoral.rep           = PDOSEoral.rat.D*BW.rat                                ## mg, amount of oral dose
  
  ex.rat.D               <- ev(ID=1, amt= DOSEoral.rep,
                               ii=tinterval.rep, addl=TDoses.rep-1,cmt="AST",replicate = FALSE)
  
  ## Exposure scenario for IV exposue to 4.2 mg/kg
  PDOSEiv.rat.A          = 4.2                                       ## Single IV dose from Luccisano et al. (2012)
  DOSEiv.rat.A           = PDOSEiv.rat.A*BW.rat                      ## amount of IV dose
  ex.iv.rat.A            <- ev(ID=1, amt= DOSEiv.rat.A, 
                           ii=tinterval.rat, addl=TDoses.rat-1, cmt="APlas_free", replicate = FALSE) # Initial dose does not bind to any proteins, so goes to the free subcompartment. 
  
  ## Exposure scenario for IV exposue to 2 mg/kg
  PDOSEiv.rat.B          = 2                                        ## mg/kg-d, single iv dose                               
  DOSEiv.rat.B           = PDOSEiv.rat.B*BW.rat                     ## mg, amount of iv dose              
  ex.iv.rat.B            <- ev(ID = 1,amt= DOSEiv.rat.B,
                            ii = tinterval.rat,addl=TDoses.rat-1,cmt="APlas_free",replicate = FALSE)
  
  ## Set up the exposure time
  tsamp.rat              = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*100,1)   ## Simulation for single oral exposure
  tsamp.rat.rep          = tgrid(0,tinterval.rat*(TDoses.rep-1)+24*100,1)   ## Siumuation 24*28 + 24*100 hours (128 days)
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  ## Exposure scenario for oral exposue to 4.2 mg/kg
  out.A <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.A, tgrid = tsamp.rat)
  
  outdf.A <- cbind.data.frame(Time   = out.A$time/24, 
                              CA     = out.A$Plasma,                             ## ug/ml,PFOS levels in plasma  
                              CL     = out.A$Liver,                              ## ug/g, PFOS levels in liver   
                              CK     = out.A$Kidney,                             ## ug/g, PFOS levels in kidney
                              Curine = (out.A$Aurine/DOSEoral.rat.A)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  out.B <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.B, tgrid = tsamp.rat)
  
  outdf.B <- cbind.data.frame(Time   = out.B$time/24, 
                              CA     = out.B$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.B$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.B$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine =(out.B$Aurine/DOSEoral.rat.B)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral exposue to 15 mg/kg
  out.C <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.C, tgrid = tsamp.rat)
  
  outdf.C <- cbind.data.frame(Time   = out.C$time/24, 
                              CA     = out.C$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.C$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.C$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine =(out.C$Aurine/DOSEoral.rat.C)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral daily exposue to 1 mg/kg-d for 28 days
  out.D <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.D, tgrid = tsamp.rat.rep)
  
  outdf.D <- cbind.data.frame(Time   = out.D$time/24, 
                              CA     = out.D$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.D$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.D$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine =(out.D$Aurine/DOSEoral.rep)*100)          ## %Dose, PFOS in urine
  
  
  ## Outiv.A: IV exposure to 4.2 mg/kg, matrix: Urine 
  outiv.A <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.iv.rat.A, tgrid=tsamp.rat)
  
  outiv.A <- cbind.data.frame(Time   = outiv.A$time/24, 
                              CA     = outiv.A$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = outiv.A$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = outiv.A$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine =(outiv.A$Aurine/DOSEiv.rat.A)*100)          ## %Dose, PFOS in urine
  
  ## Exposure scenario for IV exposue to 2 mg/kg
  outiv.B <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.iv.rat.B, tgrid=tsamp.rat)
  
  outiv.B <- cbind.data.frame(Time   = outiv.B$time/24, 
                              CA     = outiv.B$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = outiv.B$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = outiv.B$Kidney,                            ## ug/g, PFOS levels in kidney
                              Curine =(outiv.B$Aurine/DOSEiv.rat.B)*100)          ## %Dose, PFOS in urine
  
  
  
  
  if (pred) return (list( "outdf.A1"  = outdf.A, # if shows an error like this "　", just delete it.
                          "outdf.A2"  = outdf.B,
                          "outdf.A3"  = outdf.B,
                          "outdf.A4"  = outiv.A,
                          "outdf.B1"  = outdf.C, 
                          "outdf.B2"  = outdf.C,
                          "outdf.B3"  = outiv.B,
                          "outdf.B4"  = outdf.D))
  
  outdf.A1 = outdf.A[which(outdf.A$Time %in% Rat.obs$A1$Time), ]
  outdf.A2 = outdf.B[which(outdf.B$Time %in% Rat.obs$A2$Time), ]
  outdf.A3 = outdf.B[which(outdf.B$Time %in% Rat.obs$A3$Time), ]
  outdf.A4 = outiv.A[which(outiv.A$Time %in% Rat.obs$A4$Time), ]
  outdf.B1 = outdf.C[which(outdf.C$Time %in% Rat.obs$B1$Time), ]
  outdf.B2 = outdf.C[which(outdf.C$Time %in% Rat.obs$B2$Time), ]
  outdf.B3 = outiv.B[which(outdf.B$Time %in% Rat.obs$B3$Time), ]
  outdf.B4 = outdf.D[which(outdf.D$Time %in% Rat.obs$B4$Time), ]
  
  
  return (list("outdf.A1"  = outdf.A1, 
               "outdf.A2"  = outdf.A2,
               "outdf.A3"  = outdf.A3,
               "outdf.A4"  = outdf.A4,
               "outdf.B1"  = outdf.B1, 
               "outdf.B2"  = outdf.B2,
               "outdf.B3"  = outdf.B3,
               "outdf.B4"  = outdf.B4))
  
}


############################################# Monkey observed data #####################################################
## Monkey.obs:                                                                                                         # 
## A1: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: plasma                                                   #
## A2: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: urine                                                    #
## B1: Cynomolgus Monkeys oral daily dose to 0.03 mg/kg for 182 days and monitored 1 year                              #
## B2: Cynomolgus Monkeys oral daily dose to 0.15 mg/kg for 182 days and monitored 1 year                              #
## B3: Cynomolgus Monkeys oral daily dose to 0.75 mg/kg for 182 days and monitored 1 year                              #
#######################################################################################################################                                                                                                                     #


## Preidction function for monkey
pred.monkey <- function(pars,pred=FALSE) {
  
  ## Define the exposure scenario for single iv exposue to 2 mg/kg
  
  BW.monkey             = 3.5                               ## Monkey body weight
  tinterval.monkey      = 24                                ## Time interval
  TDoses.monkey         = 1                                 ## Dose times
  PDOSEiv.monkey        = 2                                 ## Single iv dose 
  DOSEiv.monkey         = PDOSEiv.monkey*BW.monkey          ## Amount of iv dose
  ex.iv.monkey          <-ev(ID = 1, amt= DOSEiv.monkey, 
                             ii = tinterval.monkey, addl = TDoses.monkey-1, cmt="APlas_free", replicate = FALSE)        
  
  ## Set up the exposure time
  tsamp.monkey.iv       = tgrid(0,tinterval.monkey*(TDoses.monkey-1)+24*200,1)   ## exposure time is 24*200 hours (200 days); 
  
  
  ## Define the exposure scenario for Repeat dose exposure scenario: 
  ## monkeys were administered 0.03, 0.15, and 0.75 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  BW.monkey.A           = 3.9                                  ## Monkey body weight
  BW.monkey.B           = 3.3                                  ## Monkey body weight
  BW.monkey.C           = 3.2                                  ## Monkey body weight
  
  tinterval.monkey      = 24                                   ## Time interval
  TDoses.monkey.rep     = 182                                  ## Dose times 182 days
  
  PDOSEoral.monkey.A    = 0.03                                 ## mg/kg; BW Oral dose
  PDOSEoral.monkey.B    = 0.15                                 ## mg/kg; BW Oral dose
  PDOSEoral.monkey.C    = 0.75                                 ## mg/kg; BW Oral dose
  
  DOSEoral.monkey.rep.A = PDOSEoral.monkey.A*BW.monkey.A       ## mg; amount of oral dose
  DOSEoral.monkey.rep.B = PDOSEoral.monkey.B*BW.monkey.B       ## mg; amount of oral dose
  DOSEoral.monkey.rep.C = PDOSEoral.monkey.C*BW.monkey.C       ## mg; amount of oral dose
  
  ex.monkey.rep.A <- ev(ID=1, amt= DOSEoral.monkey.rep.A, ii = tinterval.monkey, addl = TDoses.monkey.rep-1, cmt="AST", replicate = FALSE)
  ex.monkey.rep.B <- ev(ID=1, amt= DOSEoral.monkey.rep.B, ii = tinterval.monkey, addl = TDoses.monkey.rep-1, cmt="AST", replicate = FALSE)
  ex.monkey.rep.C <- ev(ID=1, amt= DOSEoral.monkey.rep.C, ii = tinterval.monkey, addl = TDoses.monkey.rep-1, cmt="AST", replicate = FALSE)
  
  ## Set up the exposure time
  tsamp.monkey.rep = tgrid(0,tinterval.monkey*(TDoses.monkey.rep-1)+24*565,24)  ## 182 days and simulated for 24*565 hours after dosing
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  ## Define the exposure scenario for single iv exposue to 2 mg/kg
  out.A <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.iv.monkey, tgrid = tsamp.monkey.iv)
  
  outdf.A <- cbind.data.frame(Time   = out.A$time/24, 
                              CA     = out.A$Plasma,                             ## ug/ml,PFOS levels in plasma  
                              CL     = out.A$Liver,                              ## ug/g, PFOS levels in liver   
                              CK     = out.A$Kidney,                             ## ug/g, PFOS levels in kidney
                              Curine = (out.A$Aurine/DOSEiv.monkey)*100)         ## %Dose, PFOS in urine
  
  ## Monkeys were administered 0.03 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  out.rep.A <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.monkey.rep.A, tgrid = tsamp.monkey.rep)
  
  outdf.rep.A <- cbind.data.frame(
    Time   = out.rep.A$time/24, 
    CA     = out.rep.A$Plasma,                                 ## ug/ml,PFOS levels in plasma  
    CL     = out.rep.A$Liver,                                  ## ug/g, PFOS levels in liver   
    CK     = out.rep.A$Kidney,                                 ## ug/g, PFOS levels in kidney
    Curine =(out.rep.A$Aurine/DOSEoral.monkey.rep.A)*100)      ## %Dose, PFOS in urine
  
  ## Monkeys were administered 0.15 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  out.rep.B <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.monkey.rep.B, tgrid = tsamp.monkey.rep)
  
  outdf.rep.B <- cbind.data.frame(
    Time   = out.rep.B$time/24, 
    CA     = out.rep.B$Plasma,                                 ## ug/ml,PFOS levels in plasma  
    CL     = out.rep.B$Liver,                                  ## ug/g, PFOS levels in liver   
    CK     = out.rep.B$Kidney,                                 ## ug/g, PFOS levels in kidney
    Curine =(out.rep.B$Aurine/DOSEoral.monkey.rep.B)*100)      ## %Dose, PFOS in urine
  
  ## Monkeys were administered 0.75 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  out.rep.C <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.monkey.rep.C, tgrid = tsamp.monkey.rep)
  
  outdf.rep.C <- cbind.data.frame(
    Time   = out.rep.C$time/24, 
    CA     = out.rep.C$Plasma,                                 ## ug/ml,PFOS levels in plasma  
    CL     = out.rep.C$Liver,                                  ## ug/g, PFOS levels in liver   
    CK     = out.rep.C$Kidney,                                 ## ug/g, PFOS levels in kidney
    Curine =(out.rep.C$Aurine/DOSEoral.monkey.rep.C)*100)      ## %Dose, PFOS in urine
  
  if (pred) return (list("outdf.A1"  = outdf.A, # if shows an error code "　", just delete it.
                         "outdf.A2"  = outdf.A,
                         "outdf.B1"  = outdf.rep.A,
                         "outdf.B2"  = outdf.rep.B,
                         "outdf.B3"  = outdf.rep.C))
  
  outdf.A1 = outdf.A[which(outdf.A$Time %in% Monkey.obs$A1$Time), ]
  outdf.A2 = outdf.A[which(outdf.A$Time %in% Monkey.obs$A2$Time), ]
  outdf.B1 = outdf.rep.A[which(outdf.rep.A$Time %in% Monkey.obs$B1$Time), ]
  outdf.B2 = outdf.rep.B[which(outdf.rep.B$Time %in% Monkey.obs$B2$Time), ]
  outdf.B3 = outdf.rep.C[which(outdf.rep.C$Time %in% Monkey.obs$B3$Time), ]
  
  return (list("outdf.A1"  = outdf.A1, 
               "outdf.A2"  = outdf.A2,
               "outdf.B1"  = outdf.B1, 
               "outdf.B2"  = outdf.B2,
               "outdf.B3"  = outdf.B3))
  
}

############################################# Human observed data######################################################
# Human.obs:                                                                                                          #
# Measured pooled Plasma data from 1977 - 2006 (Haugh et al., 2009)                                                   #          
#                                                                                                                     #
#######################################################################################################################

# prediction function

pred.human  <- function(pars, pred=FALSE) {
  
  # HIGH dose exposure scenario (singal dose to 2.3-3.0e-3 ug/kg/day)
  BW.human                      = 82.3
  tinterval.human               = 24
  TDoses.human                  = 365*25 # from 1976 to 2001
  
  PDOSEoral1.human              = 0.003                       # ug/kg/day Oral dose; daily dose to 2.3-3.0e-3 ug/kg/day in 1999-2000  (Loccisano et al., 2011)
  DOSEoral1.human               = PDOSEoral1.human*BW.human
  PDOSEoral2.human              = 0.0019                      # ug/kg/day Oral dose; daily dose to 1.6e-1.9e-3 ug/kg/day in 2000-2003 (Loccisano et al., 2011)
  DOSEoral2.human               = PDOSEoral2.human *BW.human
  
  # (PFOS Plasmas conc. for 1999-2000 with exposure to 2.3-3.0e-3 ug/kg/day)
  # (PFOS Plasmas conc. for 2003-2004 with exposure to 1.6-1.9e-3 ug/kg/day)
  
  oral1       <- ev(ID=1, time = 0, amt= DOSEoral1.human, ii=tinterval.human, addl = TDoses.human-1, cmt="AST", replicate = FALSE)
  oral2       <- ev(ID=1, time = 24*365*25,amt= DOSEoral2.human, ii=tinterval.human, addl = TDoses.human-1, cmt="AST", replicate = FALSE) # Initially, run high dose, then after 25 years, switch to low dose
  ex.human    <- oral1 + oral2 
  
  # set up the exposure time
  tsamp.human = tgrid(0,tinterval.human*(TDoses.human -1)+24*365*10,24*365) # Total simulation time is 25 + 10 years, sampling the data every year.
  
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  out <- 
    mod.human %>% 
    param(pars) %>%
    Req(Plasma) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d(data = ex.human, tgrid = tsamp.human)
  
  outdf <-cbind.data.frame (Time = out$time/(24*365),
                            CA   = out$Plasma)
  
  if (pred) return (outdf) # If pred=true, then run this line of code and it will generate prediction data for the entire time course. If pred=false, then run the line of the code below.
  
  outdf = outdf [which(outdf$Time %in% Human.obs$Time), ] # "%in%" overlapped the observed and predicted data at the observed time points into the same object
  
}

##################### Fig.4a Global fitting analysis ########################
## Rat
predf.rat         = pred.rat (Rat.MCMC$bestpar[1:17]) # prediction results based on best-estimated parameters from MCMC
predf.rat.A       = pred.rat (theta.Rat[1:17])

predf.rat.A1      = cbind.data.frame(
  pre.value       = predf.rat.A$outdf.A1$CA, # predf.rat.A has A1, A2, A4; outdf.A2 has CA, CL, and CK; outdf.A1 exists in different species, but you obtain this object from each specific specific, like predf.rat.A 
  obs.value       = Rat.obs$A1$CA)

predf.rat.A2      = cbind.data.frame(
  pre.value       = predf.rat.A$outdf.A2$CA,
  obs.value       = Rat.obs$A2$CA)

predf.rat.A3      = cbind.data.frame(
  pre.value       = predf.rat$outdf.A3$CA,
  obs.value       = Rat.obs$A3$CA[c(3,5:15)]) # Set the time points to match the prediction time points

predf.rat.A4    = cbind.data.frame(
  pre.value       = predf.rat.A$outdf.A4$Curine[2:14], # The first data point is an apparent outlier, so the first data point is excluded.
  obs.value       = Rat.obs$A4$Curine[2:14])

predf.rat.B1    = cbind.data.frame(
  pre.value       = predf.rat$outdf.B1$CA,
  obs.value       = Rat.obs$B1$CA)

predf.rat.B2    = cbind.data.frame(
  pre.value       = predf.rat$outdf.B2$Curine[2:11], # The first data point is an apparent outlier, so the first data point is excluded.
  obs.value       = Rat.obs$B2$Curine[2:11])

predf.rat.B3    = cbind.data.frame(
  pre.value       = predf.rat$outdf.B3$CA,
  obs.value       = Rat.obs$B3$CA)

predf.rat.B4    = cbind.data.frame(
  pre.value       = predf.rat$outdf.B4$CA,
  obs.value       = Rat.obs$B4$CA)

ep.rat          = rbind(predf.rat.A1,predf.rat.A2,predf.rat.A3,predf.rat.A4,
                        predf.rat.B1, predf.rat.B2, predf.rat.B3,predf.rat.B4) # combination of many sets of predicted vs observed data

ep.rat$Species  = rep("Rat",length(ep.rat $obs.value)) # Create a Species column, the number of rows is the same as the number of obs.value 

plot(x=log10(ep.rat$pre.value),y = log10(ep.rat$obs.value))
abline(lm(log10(ep.rat$obs.value) ~ log10(ep.rat$pre.value)))


## Mouse
predf.mouse       = pred.mouse (Mouse.MCMC$bestpar[1:17])
predf.mouse.A     = pred.mouse (theta.Mouse[1:17])

predf.mouse.A1    = cbind.data.frame(
  pre.value       = predf.mouse.A$outdf.A$CA,
  obs.value       = Mouse.obs$A1$CA)

predf.mouse.A2    = cbind.data.frame(
  pre.value       = predf.mouse.A$outdf.A$CL,
  obs.value       = Mouse.obs$A2$CL)

predf.mouse.A3    = cbind.data.frame(
  pre.value       = predf.mouse.A$outdf.A$CK,
  obs.value       = Mouse.obs$A3$CK)

predf.mouse.A4    = cbind.data.frame(
  pre.value       = predf.mouse.A$outdf.A$Curine[2:8],
  obs.value       = Mouse.obs$A4$Curine[2:8])


predf.mouse.B1    = cbind.data.frame(
  pre.value       = predf.mouse$outdf.B$CA,
  obs.value       = Mouse.obs$B1$CA)

predf.mouse.B2    = cbind.data.frame(
  pre.value       = predf.mouse$outdf.B$CL,
  obs.value       = Mouse.obs$B2$CL)

predf.mouse.B3    = cbind.data.frame(
  pre.value       = predf.mouse$outdf.B$CK,
  obs.value       = Mouse.obs$B3$CK)

predf.mouse.B4    = cbind.data.frame(
  pre.value       = predf.mouse$outdf.B$Curine [2:8], # The first data point is an apparent outlier, so the first data point is excluded.
  obs.value       = Mouse.obs$B4$Curine[2:8])

ep.mouse          = rbind(predf.mouse.A1,predf.mouse.A2,predf.mouse.A3,predf.mouse.A4,
                          predf.mouse.B1, predf.mouse.B2,predf.mouse.B3,predf.mouse.B4)

ep.mouse$Species  = rep("Mouse",length(ep.mouse$obs.value))

plot(x=log10(ep.mouse$pre.value),y = log10(ep.mouse$obs.value))
abline(lm(log10(ep.mouse$obs.value) ~ log10(ep.mouse$pre.value)))


## Monkey
predf.monkey.A    = pred.monkey (theta.Monkey[1:17])
predf.monkey      = pred.monkey (Monkey.MCMC$bestpar[1:17])

predf.monkey.A1   = cbind.data.frame(
  pre.value   = predf.monkey.A$outdf.A1$CA,
  obs.value   = Monkey.obs$A1$CA)

predf.monkey.A2   = cbind.data.frame(
  pre.value   = predf.monkey.A$outdf.A2$Curine[2:10], # The first data point is an apparent outlier, so the first data point is excluded.
  obs.value   = Monkey.obs$A2$Curine[2:10])

predf.monkey.B1   = cbind.data.frame(
  pre.value   = predf.monkey$outdf.B1$CA,
  obs.value   = Monkey.obs$B1$CA)


predf.monkey.B2   = cbind.data.frame(
  pre.value   = predf.monkey$outdf.B2$CA,
  obs.value   = Monkey.obs$B2$CA)

predf.monkey.B3   = cbind.data.frame(
  pre.value   = predf.monkey$outdf.B3$CA,
  obs.value   = Monkey.obs$B3$CA)

ep.monkey         = rbind(predf.monkey.A1,predf.monkey.A2,predf.monkey.B1,
                          predf.monkey.B2,predf.monkey.B3)

ep.monkey$Species = rep("Monkey",length(ep.monkey$obs.value))

plot(x=log(ep.monkey$pre.value),y=log(ep.monkey$obs.value))
abline(lm(log(ep.monkey$obs.value) ~ log(ep.monkey$pre.value)))

## Human
predf.human      = pred.human (theta.Human[1:17])

predf.human      = cbind.data.frame(
  pre.value      = unique(predf.human$CA),
  obs.value      = Human.obs$CA)


ep.human         = predf.human
ep.human$Species = rep("Human",length(ep.human$obs.value))

plot(x=log(ep.human$pre.value),y=log(ep.human$obs.value))
abline(lm(log(ep.human$obs.value) ~ log(ep.human$pre.value)))

############## Data of global evaluation of model fit for Figure 4a ################
## All data
ep.all =rbind(ep.mouse,ep.rat,ep.monkey,ep.human)
ep.all$log.obs = log(ep.all$obs.value,10)
ep.all$log.pre = log(ep.all$pre.value,10)
fit <- lm(log.obs ~ log.pre, data=ep.all)                   # Fit the model
ep.all$residuals = residuals(fit)                           # Save the residual values
ep.all$predicted = predict(fit)                             # Save the predicted values
ep.all$OPratio = ep.all$pre.value/ep.all$obs.value          # Estimated the predicted-to-observed ratio

p <- 
  ggplot(ep.all, aes(log.obs, log.pre)) + 
  geom_segment (aes(xend    = log.obs,                   # connect the actual data points with their corresponding predicted value    
                    yend    = predicted),
                alpha   = .2) +
  geom_point   (aes(shape   = as.factor(Species)),size = 4)  +
  #scale_colour_manual(values=c("2"="red")) +
  #scale_fill_continuous(name = "Points") +               # legend name
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black",size = 1)


p<-
  p +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x))

windowsFonts(Times=windowsFont("Times New Roman"))

p1 <- p + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="#f7f7f7"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "Observed value", y = "Predicted value")

p1

################## Figure 4b: predicted-to-observed vs. prediction plot ####################
library (DescTools)
library(ggExtra)
library(gridExtra)


p2 <-
  ggplot(ep.all, aes(predicted, log(OPratio,10))) +
  geom_hline(yintercept = log10(2),linetype = 3,color   = "red", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3,color   = "red", size =1) +
  #geom_ref_line(colour = "white", h = 0, size =0.5) +
  geom_point(color   = "#FFA488", 
             aes(shape= as.factor(Species)),size = 3) +
  geom_smooth(se = FALSE) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x))

p2 <- p2 +
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="#f7f7f7"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none')+
    labs (x = "Predicted value", y = "Predicted/Observed")

p3 <-ggMarginal(p2, type = "histogram", margins = "y",  
                yparams = list(binwidth = 0.1, fill = "#FFA488"))

grid.arrange(p1, p3, nrow = 1)

## Save the plot with high quality resolution

ggsave("Fig.4.tiff",scale = 1,
       plot = grid.arrange(p1, p3, nrow = 1),
       path = "C:/Users/weichunc/Desktop",
       width = 25, height = 12, units = "cm",dpi=320)


dev.off()

