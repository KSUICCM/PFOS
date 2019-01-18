## Load libraries
library(mrgsolve)  # Needed to run the main PBPK code
library(magrittr)  # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)   # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr)     # Needed for the pipe %>% operator
library(reshape)   # melt function to reshape the table
library(scales)    # for plotting the figure
library(gridExtra) # for plotting the figure
library(grid)      # for plotting the figure
library(lattice)   # for plotting the figure

## input PBPK model
micePBPK.code     <- readRDS (file = "micePBPK.RDS")
monkeyPBPK.code   <- readRDS (file = "monkeyPBPK.RDS")
ratPBPK.code      <- readRDS (file = "ratPBPK.RDS")
humanPBPK.code    <- readRDS (file = "humanPBPK.RDS")

## Loading PBPK model 
mod.mouse         <- mcode ("micepbpk", micePBPK.code)
mod.rat           <- mcode ("ratpbpk", ratPBPK.code)
mod.monkey        <- mcode ("monkeypbpk", monkeyPBPK.code)
mod.human         <- mcode ("humanpbpk", humanPBPK.code)


## Loading human, rat, mouse, monkey observed data
Human.obs         <- readRDS(file = "Human.obs.rds")
Rat.obs           <- readRDS(file = "Rat.obs.rds")
Mouse.obs         <- readRDS(file = "Mouse.obs.rds")
Monkey.obs        <- readRDS(file = "Monkey.obs.rds")

## Loading human, rat, mouse, monkey MCMC data (Posterior distributions)
Human.MCMC        <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC          <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC        <- readRDS(file = "mouse.MCMC.rds")
Monkey.MCMC       <- readRDS(file = "Monkey.MCMC.rds")

## Population mean (u) (Prior distributions)
theta.Human  = readRDS (file="theta.Human.Rds")
theta.Monkey = readRDS (file="theta.Monkey.Rds")
theta.Mouse  = readRDS (file="theta.Mouse.Rds")
theta.Rat    = readRDS (file="theta.Rat.Rds")

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)


## Theme
windowsFonts(Times=windowsFont("Times New Roman"))
Theme.Fig <-theme (
  #plot.background         = element_blank(),
  text                    = element_text (family = "Times",face="bold",size = 18,colour="black"),
  #panel.border            = element_rect (colour = "black", fill=NA, size=2),
  panel.background        = element_rect (fill   = "#f7f7f7"),
  #panel.grid.major.y      = element_line (size   = 0.5, colour = "grey"),
  #panel.grid.minor.y      = element_blank(),
  panel.grid.major.x      = element_blank(), 
  #panel.grid.minor.x      = element_blank(), 
  #axis.text               = element_text (size   = 18, colour = "black",face="bold"),  
  axis.title              = element_text (size   = 14, colour = "black", face = "bold"),   # label of axes
  axis.ticks.x            = element_line (size   = 0.9, colour = "black"), 
  axis.ticks.y            = element_line (size   = 0.9, colour = "black"),
  #axis.text.y             = element_blank(),
  strip.text              = element_text (size   = 12),
  legend.title            = element_text (size   = 14, face="bold"),
  legend.justification    = "none",
  legend.position         = "none",
  legend.text             = element_text (size = 14,face="bold")) 

############################################# Figure S9-S12 ###########################################################

############################################# Mouse observed data #####################################################                                                                                                                 #      
# A1: CD1 mice oral single dose of 1  mg/kg; plasma                                                                   #
# A2: CD1 mice oral single dose of 1  mg/kg; liver                                                                    #
# A3: CD1 mice oral single dose of 1  mg/kg; kidney                                                                   #
# A4: CD1 mice oral single dose of 1  mg/kg; urine                                                                    # 
# B1: CD1 mice oral single dose of 20 mg/kg; plasma                                                                   #
# B2: CD1 mice oral single dose of 20 mg/kg; liver                                                                    #
# B3: CD1 mice oral single dose of 20 mg/kg; kidney                                                                   #
# B4: CD1 mice oral single dose of 20 mg/kg; urine                                                                    #
#######################################################################################################################

## Prediction function for mouse
pred.mouse <- function(pars, pred=FALSE) {
  
  ## Define the exposure scenario for oral exposue to 1 and 20 mg/kg-d
  BW.mouse             = 0.025                                   ## Mice body weight
  tinterval.mouse      = 24                                      ## Time interval
  TDoses.mouse         = 1                                       ## Dose times
  PDOSEoral.mouse.A    = 1                                       ## 1  mg/kg-d, Single oral dose
  PDOSEoral.mouse.B    = 20                                      ## 20 mg/kg-d, Single oral dose
  
  DOSEoral.mouse.A     = PDOSEoral.mouse.A*BW.mouse              ## amount of oral dose
  DOSEoral.mouse.B     = PDOSEoral.mouse.B*BW.mouse              ## amount of oral dose
  
  ex.mouse.A           <- ev (ID = 1, amt= DOSEoral.mouse.A, 
                              ii = tinterval.mouse, addl=TDoses.mouse-1, cmt="AST", replicate = FALSE)
  ex.mouse.B           <- ev (ID = 1, amt= DOSEoral.mouse.B, 
                              ii = tinterval.mouse, addl=TDoses.mouse-1, cmt="AST", replicate = FALSE)
  
  
  ## set up the exposure time
  tsamp.mouse          = tgrid(0,tinterval.mouse*(TDoses.mouse-1)+24*200,1)     ## Simulation for 24*200 hours (200 days)
  
  
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  ## Get a prediction
  out.A <- 
    mod.mouse %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.mouse.A, tgrid = tsamp.mouse)
  
  outdf.A <- cbind.data.frame(Time   = out.A$time/24, 
                              CA     = out.A$Plasma,                             ## ug/ml,PFOS levels in plasma  
                              CL     = out.A$Liver,                              ## ug/g, PFOS levels in liver   
                              CK     = out.A$Kidney,                             ## ug/g, PFOS in kidney
                              Curine = (out.A$Aurine/DOSEoral.mouse.A)*100)      ## %Dose, PFOS in urine
  out.B <- 
    mod.mouse %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d (data = ex.mouse.B, tgrid = tsamp.mouse)
  
  outdf.B <- cbind.data.frame(Time   = out.B$time/24, 
                              CA     = out.B$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.B$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.B$Kidney,                            ## ug/g, PFOS in kidney
                              Curine = (out.B$Aurine/DOSEoral.mouse.B)*100)     ## %Dose, PFOS in urine
  
  if (pred) return (list("outdf.A"= outdf.A, "outdf.B" = outdf.B))
  
  
  
   
  outdf.A = outdf.A[which(outdf.A$Time %in% Mouse.obs$A1$Time),] # Which indicates to extract data when predicted values at the time points when there were observed data.
  outdf.B = outdf.B[which(outdf.B$Time %in% Mouse.obs$B1$Time),]
  
  return (list("outdf.A"= outdf.A, "outdf.B" = outdf.B))
  
}


############################################# Rat data ################################################################
# Rat.obs:                                                                                                            #    
# A1: SD rat oral single dose to 4.2 mg/kg, matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #                  
# A2: SD rat oral single dose to 2   mg/kg, matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #
# A3: SD rat oral single dose to 2   mg/kg, matrix: plasma, data from kim et al. (2016)                               #
# A4: SD rat iv exposure to 4.2 mg/kg,      matrix: urine, data from Johnson et al. (1979), Luccisano et al. (2012);  #
# B1: SD rat oral single dose to 15 mg/kg,  matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #
# B2: SD rat oral single dose to 15 mg/kg,  matrix: urine, data from Chang et al. (2012)                              #
# B3: SD rat iv single dose to 2 mg/kg,     matrix: plasma, data from kim et al. (2016)                               #
# B4: SD rat oral daily dose to 1 mg/kg for 28 days,matrix: plasma,data from Change et al. (2012)                     #
#######################################################################################################################


## Preiction function
pred.rat <- function(pars,pred=FALSE) {
  
  ## Exposure scenario for oral exposue to 4.2 mg/kg
  
  BW.rat                 = 0.3                                       # Rat body weight
  tinterval.rat          = 24                                        # Time interval
  TDoses.rat             = 1                                         # Dose times
  PDOSEoral.rat.A        = 4.2                                       # Single oral dose from Luccisano et al. (2012)
  DOSEoral.rat.A         = PDOSEoral.rat.A*BW.rat                    # amount of oral dose
  ex.rat.A               <- ev (ID = 1, amt= DOSEoral.rat.A, 
                                ii = tinterval.rat, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  PDOSEoral.rat.B        = 2                                         ## Single oral dose from Luccisano et al. (2012)
  DOSEoral.rat.B         = PDOSEoral.rat.B*BW.rat                    ## amount of oral dose
  ex.rat.B               <- ev (ID=1, amt= DOSEoral.rat.B, 
                                ii=tinterval.rat, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  
  ## Exposure scenario for oral exposue to 15 mg/kg
  PDOSEoral.rat.C        = 15                                        ## mg/kg-d, Single oral dose                           
  DOSEoral.rat.C         = PDOSEoral.rat.C*BW.rat                    ## mg, amount of oral dose                   
  ex.rat.C               <- ev(ID=1,amt= DOSEoral.rat.C, 
                               ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  
  ## Exposure scenario for oral daily exposue to 1 mg/kg-d for 28 days
  tinterval.rep          = 24                                                    ## hr,    time interval
  TDoses.rep             = 28                                                    ## times, Dose times
  
  PDOSEoral.rat.D        = 1                                                     ## mg/kg, BW Oral dose
  DOSEoral.rep           = PDOSEoral.rat.D*BW.rat                                ## mg, amount of oral dose
  
  ex.rat.D               <- ev(ID=1, amt= DOSEoral.rep,
                               ii=tinterval.rep, addl=TDoses.rep-1,cmt="AST",replicate = FALSE)
  
  ## Exposure scenario for IV exposue to 4.2 mg/kg
  PDOSEiv.rat.A      = 4.2                                       ## Single IV dose from Johnson et al. (1979)
  DOSEiv.rat.A       = PDOSEiv.rat.A*BW.rat                      ## amount of IV dose
  ex.iv.rat.A        <- ev(ID=1, amt= DOSEiv.rat.A, 
                           ii=tinterval.rat, addl=TDoses.rat-1, cmt="APlas_free", replicate = FALSE)
  
  ## Exposure scenario for IV exposue to 2 mg/kg
  PDOSEiv.rat.B       = 2                                        ## mg/kg-d, Single iv dose                               
  DOSEiv.rat.B        = PDOSEiv.rat.B*BW.rat                     ## mg, amount of iv dose              
  ex.iv.rat.B         <- ev(ID = 1,amt= DOSEiv.rat.B,
                            ii = tinterval.rat,addl=TDoses.rat-1,cmt="APlas_free",replicate = FALSE)
  
  
  ## set up the exposure time
  tsamp.rat           = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*100,1)
  tsamp.rat.rep       = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*100,1)   ## Siumuation 24*28 + 24*100 hours (128 days)
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  ## Exposure scenario for oral exposue to 4.2 mg/kg
  out.A <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.A, tgrid = tsamp.rat)
  
  outdf.A <- cbind.data.frame(Time   = out.A$time/24, 
                              CA     = out.A$Plasma,                             ## ug/ml,PFOS levels in plasma  
                              CL     = out.A$Liver,                              ## ug/g, PFOS levels in liver   
                              CK     = out.A$Kidney,                             ## ug/g, PFOS in kidney
                              Curine = (out.A$Aurine/DOSEoral.rat.A)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  out.B <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.B, tgrid = tsamp.rat)
  
  outdf.B <- cbind.data.frame(Time   = out.B$time/24, 
                              CA     = out.B$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.B$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.B$Kidney,                            ## ug/g, PFOS in kidney
                              Curine =(out.B$Aurine/DOSEoral.rat.B)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral exposue to 15 mg/kg
  out.C <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.C, tgrid = tsamp.rat)
  
  outdf.C <- cbind.data.frame(Time   = out.C$time/24, 
                              CA     = out.C$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.C$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.C$Kidney,                            ## ug/g, PFOS in kidney
                              Curine =(out.C$Aurine/DOSEoral.rat.C)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral daily exposue to 1 mg/kg-d for 28 days
  out.D <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.rat.D, tgrid = tsamp.rat.rep)
  
  outdf.D <- cbind.data.frame(Time   = out.D$time/24, 
                              CA     = out.D$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = out.D$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = out.D$Kidney,                            ## ug/g, PFOS in kidney
                              Curine =(out.D$Aurine/DOSEoral.rep)*100)        ## %Dose, PFOS in urine
  
  
  ## Outiv.A: IV exposure to 4.2 mg/kg, matrix: Plasma, urine 
  outiv.A <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.iv.rat.A, tgrid=tsamp.rat)
  
  outiv.A <- cbind.data.frame(Time   = outiv.A$time/24, 
                              CA     = outiv.A$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = outiv.A$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = outiv.A$Kidney,                            ## ug/g, PFOS in kidney
                              Curine =(outiv.A$Aurine/DOSEiv.rat.A)*100)          ## %Dose, PFOS in urine
  
  ## Exposure scenario for IV exposue to 2 mg/kg
  outiv.B <- 
    mod.rat %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.iv.rat.B, tgrid=tsamp.rat)
  
  outiv.B <- cbind.data.frame(Time   = outiv.B$time/24, 
                              CA     = outiv.B$Plasma,                            ## ug/ml,PFOS levels in plasma  
                              CL     = outiv.B$Liver,                             ## ug/g, PFOS levels in liver   
                              CK     = outiv.B$Kidney,                            ## ug/g, PFOS in kidney
                              Curine =(outiv.B$Aurine/DOSEiv.rat.B)*100)          ## %Dose, PFOS in urine
  
  
  
  
  if (pred) return (list("outdf.A1"  = outdf.A, 
                          "outdf.A2"  = outdf.B,
                          "outdf.A3"  = outdf.B,
                          "outdf.A4"  = outiv.A,
                          "outdf.B1"  = outdf.C, 
                          "outdf.B2"  = outdf.C,
                          "outdf.B3"  = outiv.B,
                          "outdf.B4"  = outdf.D))

  outdf.A1 = outdf.A[which(outdf.A$Time %in% Rat.obs$A1$Time), ] # Which indicates to extract data when predicted values at the time points when there were observed data.
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


############################################# Monkey data ################################################################
## Monkey.obs:                                                                                                           # 
## A1: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: plasma; chang et al., 2010                                 #
## A2: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: urine; chang et al., 2010                                  #
## B1: Cynomolgus Monkeys oral daily dose to 0.03 mg/kg for 182 days and monitored 1 year; Secat et al., 2002            #
## B2: Cynomolgus Monkeys oral daily dose to 0.15 mg/kg for 182 days and monitored 1 year; Secat et al., 2002            #
## B3: Cynomolgus Monkeys oral daily dose to 0.75 mg/kg for 182 days and monitored 1 year; Secat et al., 2002            #
##########################################################################################################################                                                                                                                     #

## Preidction function for monkey
pred.monkey <- function(pars,pred=FALSE) {
  ## Define the exposure scenario for single iv exposue to 2 mg/kg
  
  BW.monkey             = 3.5                               ## Monkey body weight
  tinterval.monkey      = 24                                ## Time interval
  TDoses.monkey         = 1                                 ## Dose times
  PDOSEiv.monkey        = 2                                 ## Single iv dose 
  DOSEiv.monkey         = PDOSEiv.monkey*BW.monkey          ## amount of iv dose
  ex.iv.monkey          <-ev(ID = 1, amt= DOSEiv.monkey, 
                             ii = tinterval.monkey, addl = TDoses.monkey-1, cmt="APlas_free", replicate = FALSE)        
  
  ## set up the simulation time
  tsamp.monkey.iv       = tgrid(0,tinterval.monkey*(TDoses.monkey-1)+24*200,1)   ## simulation time is 24*200 hours (200 days); 
  
  
  ## Define the exposure scenario for Repeat dose Exposure scenario: 
  ## monkeys were administered 0.03, 0.15, and 0.75 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  BW.monkey.A           = 3.9                                  ## Monkey body weight
  BW.monkey.B           = 3.3
  BW.monkey.C           = 3.2
  
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
  
  ## set up the exposure time
  tsamp.monkey.rep = tgrid(0,tinterval.monkey*(TDoses.monkey.rep-1)+24*565,24)  ## 182 days and simulated for 24*565 hours after dosing
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  ## Define the exposure scenario for single iv exposue to 2 mg/kg-d
  out.A <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.iv.monkey, tgrid = tsamp.monkey.iv)
  
  outdf.A <- cbind.data.frame(Time   = out.A$time/24, 
                              CA     = out.A$Plasma,                             ## ug/ml,PFOS levels in plasma  
                              CL     = out.A$Liver,                              ## ug/g, PFOS levels in liver   
                              CK     = out.A$Kidney,                             ## ug/g, PFOS in kidney
                              Curine = (out.A$Aurine/DOSEiv.monkey)*100)         ## %Dose, PFOS in urine
  
  ## Monkeys were administered 0.03 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  out.rep.A <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.monkey.rep.A, tgrid = tsamp.monkey.rep)
  
  outdf.rep.A <- cbind.data.frame(
                                Time   = out.rep.A$time/24, 
                                CA     = out.rep.A$Plasma,                                 ## ug/ml,PFOS levels in plasma  
                                CL     = out.rep.A$Liver,                                  ## ug/g, PFOS levels in liver   
                                CK     = out.rep.A$Kidney,                                 ## ug/g, PFOS in kidney
                                Curine =(out.rep.A$Aurine/DOSEoral.monkey.rep.A)*100)      ## %Dose, PFOS in urine
  
  ## Monkeys were administered 0.15 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  out.rep.B <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.monkey.rep.B, tgrid = tsamp.monkey.rep)
  
  outdf.rep.B <- cbind.data.frame(
                                 Time   = out.rep.B$time/24, 
                                 CA     = out.rep.B$Plasma,                                 ## ug/ml,PFOS levels in plasma  
                                 CL     = out.rep.B$Liver,                                  ## ug/g, PFOS levels in liver   
                                 CK     = out.rep.B$Kidney,                                 ## ug/g, PFOS in kidney
                                 Curine =(out.rep.B$Aurine/DOSEoral.monkey.rep.B)*100)      ## %Dose, PFOS in urine
  
  ## Monkeys were administered 0.75 mg/kg PFOS orally for 182 days and monitored for a year after dosing
  out.rep.C <- 
    mod.monkey  %>% 
    param(pars) %>%
    Req (Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d (data = ex.monkey.rep.C, tgrid = tsamp.monkey.rep)
  
  outdf.rep.C <- cbind.data.frame(
                                Time   = out.rep.C$time/24, 
                                CA     = out.rep.C$Plasma,                                 ## ug/ml,PFOS levels in plasma  
                                CL     = out.rep.C$Liver,                                  ## ug/g, PFOS levels in liver   
                                CK     = out.rep.C$Kidney,                                 ## ug/g, PFOS in kidney
                                Curine =(out.rep.C$Aurine/DOSEoral.monkey.rep.C)*100)      ## %Dose, PFOS in urine
  
  
  
  if (pred) return (list( "outdf.A1"  = outdf.A, 
                          "outdf.A2"  = outdf.A,
                          "outdf.B1"  = outdf.rep.A,
                          "outdf.B2"  = outdf.rep.B,
                          "outdf.B3"  = outdf.rep.C))
                

  outdf.A1 = outdf.A[which(outdf.A$Time %in% Monkey.obs$A1$Time), ] # Which indicates to extract data when predicted values at the time points when there were observed data.
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



############################################# Human data ###########################################################
# Human.obs:                                                                                                          #
# Measured pooled Plasma from 1977 - 2006 (Haugh et al., 2009)                                                        #          
# This is the only human study used in the calibration and optimization. Other datasets were used in the evaluation.                                                                                                                   #
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
  tsamp.human = tgrid(0,tinterval.human*(TDoses.human -1)+25*365*10,24*365) # Total simulation time is 25 + 10 years, sampling the data every year.
  
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  out <- 
    mod.human %>% 
    param(pars) %>%
    Req(Plasma) %>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d(data = ex.human, tgrid = tsamp.human)
  
  outdf <-cbind.data.frame (Time = out$time/(24*365),
                            CA   = out$Plasma)
  
  if (pred) return (outdf) # If pred=true, then run this line of code and it will generate prediction data for the entire time course. If pred=false, then run the line of the code below.
  
  outdf = outdf [which(outdf$Time %in% Human.obs$Time), ] # "%in%" overlapped the observed and predicted data at the observed time points into the same object
  
}


###############################  Time course plot #################################

## Mouse: Chang et al. (2012)
############################################# Mouse observed data #####################################################                                                                                                                 #      
# A1: CD1 mice oral single dose of 1  mg/kg; plasma                                                                   #
# A2: CD1 mice oral single dose of 1  mg/kg; liver                                                                    #
# A3: CD1 mice oral single dose of 1  mg/kg; kidney                                                                   #
# A4: CD1 mice oral single dose of 1  mg/kg; urine                                                                    # 
# B1: CD1 mice oral single dose of 20 mg/kg; plasma                                                                   #
# B2: CD1 mice oral single dose of 20 mg/kg; liver                                                                    #
# B3: CD1 mice oral single dose of 20 mg/kg; kidney                                                                   #
# B4: CD1 mice oral single dose of 20 mg/kg; urine                                                                    #
#######################################################################################################################

## Simulation of dose group A (1 mg/kg/day) and group B (20 mg/kg/day)
## Make a table frame
Sim.mouse.A = pred.mouse (theta.Mouse[1:17],pred=TRUE)$outdf.A # Calibration results
Sim.mouse.B = pred.mouse (Mouse.MCMC$bestpar[1:17],pred=TRUE)$outdf.B # Optimization results


mouse.sim.A =
  cbind.data.frame (Time     = Sim.mouse.A$Time, 
                    Plasma   = Sim.mouse.A$CA, 
                    Liver    = Sim.mouse.A$CL,
                    Kidney   = Sim.mouse.A$CK,
                    Urine    = Sim.mouse.A$Curine)

mouse.sim.B =
  cbind.data.frame (Time     = Sim.mouse.B$Time, 
                    Plasma   = Sim.mouse.B$CA, 
                    Liver    = Sim.mouse.B$CL,
                    Kidney   = Sim.mouse.B$CK,
                    Urine    = Sim.mouse.B$Curine)

mouse.sim.A1               = melt(mouse.sim.A[,2:5]) # The format of the original mouse.sim.A is Time, CA, CL, CK, Curine. After melt function, the format is Tissue, Value
mouse.sim.A1$Time          = rep(Sim.mouse.A$Time,4) # Add the column of Time as the third column
names(mouse.sim.A1)        = c("Tissue","value","Time") # Rename the columns of the new data format

mouse.sim.B1               = melt(mouse.sim.B[,2:5])
mouse.sim.B1$Time          = rep(Sim.mouse.B$Time,4)
names(mouse.sim.B1)        = c("Tissue","value","Time")

mouse.obs.A = 
  cbind.data.frame  (Time    = Mouse.obs$A1$Time,
                     Plasma  = Mouse.obs$A1$CA,
                     Liver   = Mouse.obs$A2$CL,
                     Kidney  = Mouse.obs$A3$CK,
                     Urine   = Mouse.obs$A4$Curine)

mouse.obs.A1               = melt (mouse.obs.A[,2:5])
mouse.obs.A1$Time          = rep  (mouse.obs.A$Time,4)
names(mouse.obs.A1)        = c("Tissue","value","Time")

mouse.obs.B = 
  cbind.data.frame (Time     = Mouse.obs$B1$Time,
                    Plasma   = Mouse.obs$B1$CA,
                    Liver    = Mouse.obs$B2$CL,
                    Kidney   = Mouse.obs$B3$CK,
                    Urine    = Mouse.obs$B4$Curine)

mouse.obs.B1               = melt (mouse.obs.B[,2:5])
mouse.obs.B1$Time          = rep (mouse.obs.B$Time,4)
names(mouse.obs.B1)        = c("Tissue","value","Time")


## plot 
plot.mouse.A =
  ggplot() +
  #geom_line  (data = mouse.sim.B1, aes(x = Time, y = value), color = "blue",size=1.5) +
  #geom_point (data = mouse.obs.B1, aes(x = Time, y = value), color = "blue") +
  geom_line  (data = mouse.sim.A1, aes(x = Time, y = value, group = Tissue),color = "red",size=1.2) +
  geom_point (data = mouse.obs.A1, aes(x = Time, y = value,shape = Tissue), color = "black",size=2.5)+
  facet_wrap ( ~ Tissue, ncol=1,scale="free")

# plot.mouse.A =  plot.mouse.A + 
#   annotation_logticks(sides = "l",size=1.2) + 
#   scale_y_log10(breaks=c(0,10,100),expand = c(0,0),limits = c(1e-1,1e3)) + xlim(c(0,150))
# 

plot.mouse.A = plot.mouse.A + Theme.Fig + labs(x="",y="")
plot.mouse.A


plot.mouse.B =
  ggplot() +
  geom_line  (data = mouse.sim.B1, aes(x = Time, y = value, group=Tissue),color = "red",size=1.2) +
  geom_point (data = mouse.obs.B1, aes(x = Time, y = value,shape=Tissue), color = "black",size=2.5) +
  #geom_line  (data = mouse.sim.A1, aes(x = Time, y = value), color = "red") +
  #geom_point (data = mouse.obs.A1, aes(x = Time, y = value), color = "red") +
  facet_wrap ( ~ Tissue, ncol=1,scale="free")

# plot.mouse.B =  plot.mouse.B + 
#   annotation_logticks(sides = "l",size=1.2) + 
#   scale_y_log10(breaks=c(0,10,100),expand = c(0,0),limits = c(1e-1,1e2)) + xlim(c(0,150))


plot.mouse.B = plot.mouse.B + Theme.Fig + labs(x="",y="")
plot.mouse.B

grid.arrange(plot.mouse.A, plot.mouse.B, ncol = 2)

### Save figure #####
ggsave("Fig.S9.tiff",scale = 1,
       plot = grid.arrange(plot.mouse.A, plot.mouse.B, ncol = 2),
       path = "C:/Users/weichunc/Desktop",
       width = 20, height = 25, units = "cm",dpi=320)


############################################ RAT plot #################################################################
# A1: SD rat oral single dose to 4.2 mg/kg, matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #                   
# A2: SD rat oral single dose to 2   mg/kg, matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #
# A3: SD rat oral single dose to 2   mg/kg, matrix: plasma, data from kim et al. (2016)                               #
# A4: SD rat iv exposure to 4.2 mg/kg,      matrix: urine, data from Johnson et al. (1979), Luccisano et al. (2012);  #
# B1: SD rat oral single dose to 15 mg/kg,  matrix: plasma, data from Luccisano et al. (2012); 3M unpublished data    #
# B2: SD rat oral single dose to 15 mg/kg,  matrix: urine, data from Chang et al. (2012)                              #
# B3: SD rat iv single dose to 2 mg/kg,     matrix: plasma, data from kim et al. (2016)                               #
# B4: SD rat oral daily dose to 1 mg/kg for 29 days,matrix: plasma,data from Change et al. (2012)                     #
#######################################################################################################################

## Generate the data needed for the figure
Sim.rat.A1  = pred.rat (theta.Rat[1:17],pred=TRUE)$outdf.A1
Sim.rat.A2  = pred.rat (theta.Rat[1:17],pred=TRUE)$outdf.A2
Sim.rat.A3  = pred.rat (theta.Rat[1:17],pred=TRUE)$outdf.A3
Sim.rat.A4  = pred.rat (theta.Rat[1:17],pred=TRUE)$outdf.A4
Sim.rat.B1  = pred.rat (Rat.MCMC$bestpar[1:17],pred=TRUE)$outdf.B1
Sim.rat.B2  = pred.rat (Rat.MCMC$bestpar[1:17],pred=TRUE)$outdf.B2
Sim.rat.B3  = pred.rat (Rat.MCMC$bestpar[1:17],pred=TRUE)$outdf.B3
Sim.rat.B4  = pred.rat (Rat.MCMC$bestpar[1:17],pred=TRUE)$outdf.B4

A1= c("Oral single dose: 4.2 (mg/kg-day)/Plasma/3M unpublished data")
A2= c("Oral single dose: 2 (mg/kg-day)/Plasma/ 3M unpublished data")
A3= c("Oral single dose: 2 (mg/kg-day)/Plasma/ Kim et al. (2016)")
A4= c("IV dose: 4.2 (mg/kg-day)/ Urine/ Johnson et al. (1979)")
B1= c("Oral single dose: 15 (mg/kg-day)/ Plasma/ 3M unpublished data")
B2= c("Oral single dose: 15 (mg/kg-day)/ Urine/ Chang et al. (2012)")
B3= c("IV dose: 2 (mg/kg-day)/ Plasma/ Kim et al. (2016)")
B4= c("Oral daily dose: 1 mg/kg-day/ Plasma/ Chang et al. (2012)")



rat.sim.A1  = melt(cbind.data.frame (Data = rep(A1,length(Sim.rat.A1$Time)), CA = Sim.rat.A1$CA))
rat.sim.A2  = melt(cbind.data.frame (Data = rep(A2,length(Sim.rat.A2$Time)), CA = Sim.rat.A2$CA)) 
rat.sim.A3  = melt(cbind.data.frame (Data = rep(A3,length(Sim.rat.A3$Time)), CA = Sim.rat.A3$CA))
rat.sim.A4  = melt(cbind.data.frame (Data = rep(A4,length(Sim.rat.A4$Time)), Curine = Sim.rat.A4$Curine)) 
rat.sim.B1  = melt(cbind.data.frame (Data = rep(B1,length(Sim.rat.B1$Time)), CA = Sim.rat.B1$CA))
rat.sim.B2  = melt(cbind.data.frame (Data = rep(B2,length(Sim.rat.B2$Time)), Curine = Sim.rat.B2$Curine))
rat.sim.B3  = melt(cbind.data.frame (Data = rep(B3,length(Sim.rat.B3$Time)), CA = Sim.rat.B3$CA)) 
rat.sim.B4  = melt(cbind.data.frame (Data = rep(B4,length(Sim.rat.B4$Time)), CA = Sim.rat.B4$CA))

rat.sim.A1$Time = Sim.rat.A1$Time 
rat.sim.A2$Time = Sim.rat.A2$Time 
rat.sim.A3$Time = Sim.rat.A3$Time 
rat.sim.A4$Time = Sim.rat.A4$Time 
rat.sim.B1$Time = Sim.rat.B1$Time 
rat.sim.B2$Time = Sim.rat.B2$Time 
rat.sim.B3$Time = Sim.rat.B3$Time 
rat.sim.B4$Time = Sim.rat.B4$Time 


rat.obs.A1 = 
  melt(cbind.data.frame   (CA     = Rat.obs$A1$CA,
                           Data   = rep(A1,length(Rat.obs$A1$Time))))

rat.obs.A2 = 
  melt(cbind.data.frame   (
    CA     = Rat.obs$A2$CA,
    Data   = rep(A2,length(Rat.obs$A2$Time))))

rat.obs.A3 = 
  melt(cbind.data.frame   (
    CA     = Rat.obs$A3$CA,
    Data   = rep(A3,length(Rat.obs$A3$Time))))

rat.obs.A4 = 
  melt(cbind.data.frame   (
    Curine     = Rat.obs$A4$Curine,
    Data   = rep(A4,length(Rat.obs$A4$Time))))


rat.obs.B1 = 
  melt(cbind.data.frame   (
    CA     = Rat.obs$B1$CA,
    Data   = rep(B1,length(Rat.obs$B1$Time))))

rat.obs.B2 = 
  melt(cbind.data.frame   (
    Curine     = Rat.obs$B2$Curine,
    Data   = rep(B2,length(Rat.obs$B2$Time))))

rat.obs.B3 = 
  melt(cbind.data.frame   (
    CA     = Rat.obs$B3$CA,
    Data   = rep(B3,length(Rat.obs$B3$Time))))

rat.obs.B4 = 
  melt(cbind.data.frame   (
    CA     = Rat.obs$B4$CA,
    Data   = rep(B4,length(Rat.obs$B4$Time))))


rat.obs.A1$Time = Rat.obs$A1$Time 
rat.obs.A2$Time = Rat.obs$A2$Time 
rat.obs.A3$Time = Rat.obs$A3$Time 
rat.obs.A4$Time = Rat.obs$A4$Time 
rat.obs.B1$Time = Rat.obs$B1$Time 
rat.obs.B2$Time = Rat.obs$B2$Time 
rat.obs.B3$Time = Rat.obs$B3$Time 
rat.obs.B4$Time = Rat.obs$B4$Time 

rat.sim.A = rbind(rat.sim.A1,rat.sim.A2,rat.sim.A3,rat.sim.A4)

rat.sim.B = rbind(rat.sim.B1,rat.sim.B2,rat.sim.B3,rat.sim.B4)

rat.obs.A = rbind(rat.obs.A1,rat.obs.A2,rat.obs.A3,rat.obs.A4)

rat.obs.B = rbind(rat.obs.B1,rat.obs.B2,rat.obs.B3,rat.obs.B4)


## plot
plot.rat.A =
  ggplot() +
  geom_line  (data = rat.sim.A, aes(x = Time, y = value,group = Data), color = "red", size=1.2) +
  geom_point (data = rat.obs.A, aes(x = Time, y = value,group = Data), color = "black",size=2.5) +
  facet_wrap ( ~ Data, ncol=1,scales="free")

# plot.rat.A =  plot.rat.A + 
#   annotation_logticks(sides = "l") +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = c(1e-1,1e2))+ xlim(c(0,100)) 

plot.rat.A = plot.rat.A + Theme.Fig + labs(x="",y="")


plot.rat.B =
  ggplot() +
  geom_line  (data = rat.sim.B, aes(x = Time, y = value,group = Data), color = "red", size=1.2) +
  geom_point (data = rat.obs.B, aes(x = Time, y = value,group = Data), color = "black",size=2.5) +
  facet_wrap ( ~ Data, ncol=1,scales="free")

# plot.rat.B =  plot.rat.B + 
#   annotation_logticks(sides = "l") +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = c(1e-2,1e2))+ xlim(c(0,100)) 

plot.rat.B = plot.rat.B + Theme.Fig + labs(x="",y="")
grid.arrange(plot.rat.A, plot.rat.B, ncol = 2)

### Save figure #####
ggsave("Fig.S10.tiff",scale = 1,
       plot = grid.arrange(plot.rat.A, plot.rat.B, ncol = 2),
       path = "C:/Users/weichunc/Desktop",
       width = 30, height = 25, units = "cm",dpi=320)


############################################# Monkey ##############################################
## A1: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: plasma                              #
## A2: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: urine                               #
## B1: Cynomolgus Monkeys oral daily dose to 0.03 mg/kg for 182 days and monitored 1 year         #
## B2: Cynomolgus Monkeys oral daily dose to 0.15 mg/kg for 182 days and monitored 1 year         #
## B3: Cynomolgus Monkeys oral daily dose to 0.75 mg/kg for 182 days and monitored 1 year         #
###################################################################################################

## Generate the data needed for the figure
Sim.monkey.A1  = pred.monkey (theta.Monkey[1:17],pred=TRUE)$outdf.A1
Sim.monkey.A2  = pred.monkey (theta.Monkey[1:17],pred=TRUE)$outdf.A2
Sim.monkey.B1  = pred.monkey (Monkey.MCMC$bestpar[1:17],pred=TRUE)$outdf.B1
Sim.monkey.B2  = pred.monkey (Monkey.MCMC$bestpar[1:17],pred=TRUE)$outdf.B2
Sim.monkey.B3  = pred.monkey (Monkey.MCMC$bestpar[1:17],pred=TRUE)$outdf.B3

A1= c("IV single dose: 2 (mg/kg-day)/Plasma/ Chang et al. (2012)")
A2= c("IV single dose: 2 (mg/kg-day)/Urine/ Chang et al. (2012)")
B1= c("Oral daily dose: 0.03 (mg/kg-day)/Plasma/ Seacat et al. (2002)")
B2= c("Oral daily dose: 0.15 (mg/kg-day)/Plasma/ Seacat et al. (2002)")
B3= c("Oral daily dose: 0.75 (mg/kg-day)/Plasma/ Seacat et al. (2002)")

monkey.sim.A1  = melt(cbind.data.frame (Data = rep(A1,length(Sim.monkey.A1$Time)), CA         = Sim.monkey.A1$CA))
monkey.sim.A2  = melt(cbind.data.frame (Data = rep(A2,length(Sim.monkey.A2$Time)), Curine     = Sim.monkey.A2$Curine)) 
monkey.sim.B1  = melt(cbind.data.frame (Data = rep(B1,length(Sim.monkey.B1$Time)), CA         = Sim.monkey.B1$CA))
monkey.sim.B2  = melt(cbind.data.frame (Data = rep(B2,length(Sim.monkey.B2$Time)), CA         = Sim.monkey.B2$CA))
monkey.sim.B3  = melt(cbind.data.frame (Data = rep(B3,length(Sim.monkey.B3$Time)), CA         = Sim.monkey.B3$CA)) 

monkey.sim.A1$Time = Sim.monkey.A1$Time 
monkey.sim.A2$Time = Sim.monkey.A2$Time 
monkey.sim.B1$Time = Sim.monkey.B1$Time 
monkey.sim.B2$Time = Sim.monkey.B2$Time 
monkey.sim.B3$Time = Sim.monkey.B3$Time 


monkey.obs.A1 = 
  melt(cbind.data.frame   (
    CA     = Monkey.obs$A1$CA,
    Data   = rep(A1,length(Monkey.obs$A1$Time))))

monkey.obs.A2 = 
  melt(cbind.data.frame   (
    Curine = Monkey.obs$A2$Curine,
    Data   = rep(A2,length(Monkey.obs$A2$Time))))

monkey.obs.B1 = 
  melt(cbind.data.frame   (
    CA     = Monkey.obs$B1$CA,
    Data   = rep(B1,length(Monkey.obs$B1$Time))))

monkey.obs.B2 = 
  melt(cbind.data.frame   (
    CA     = Monkey.obs$B2$CA,
    Data   = rep(B2,length(Monkey.obs$B2$Time))))

monkey.obs.B3 = 
  melt(cbind.data.frame   (
    CA     = Monkey.obs$B3$CA,
    Data   = rep(B3,length(Monkey.obs$B3$Time))))


monkey.obs.A1$Time = Monkey.obs$A1$Time 
monkey.obs.A2$Time = Monkey.obs$A2$Time 
monkey.obs.B1$Time = Monkey.obs$B1$Time 
monkey.obs.B2$Time = Monkey.obs$B2$Time 
monkey.obs.B3$Time = Monkey.obs$B3$Time 


monkey.sim.A = rbind(monkey.sim.A1,monkey.sim.A2)
monkey.sim.B = rbind(monkey.sim.B1,monkey.sim.B2,monkey.sim.B3)

monkey.obs.A = rbind(monkey.obs.A1,monkey.obs.A2)
monkey.obs.B = rbind(monkey.obs.B1,monkey.obs.B2,monkey.obs.B3)


plot.monkey.A =
  ggplot() +
  geom_line  (data = monkey.sim.A, aes(x = Time, y = value, group = variable), color = "red", size=1.2) +
  geom_point (data = monkey.obs.A, aes(x = Time, y = value, group = variable, shape = variable), color = "black",size=2.5)+ 
  facet_wrap ( ~ Data, ncol = 1,scales="free")

# plot.monkey.A =  plot.monkey.A + 
#   annotation_logticks(sides = "l") +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = c(1e-1,1e3))+ xlim(c(0,150))

plot.monkey.A = plot.monkey.A + Theme.Fig + labs(x="",y="")
plot.monkey.A 

plot.monkey.B =
  ggplot() +
  geom_line  (data = monkey.sim.B, aes(x = Time, y = value, group = Data), color = "red", size=1.2) +
  geom_point (data = monkey.obs.B, aes(x = Time, y = value, group = Data, shape = Data), color = "black",size=2.5)+ 
  facet_wrap ( ~ Data, ncol = 1,scales="free")

# plot.monkey.B =  plot.monkey.B + 
#   annotation_logticks(sides = "l") +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = c(1e-1,1e3))+ xlim(c(0,600))

plot.monkey.B = plot.monkey.B + Theme.Fig + labs(x="",y="")

grid.arrange(plot.monkey.A, plot.monkey.B, ncol = 2)


### Save figure #####
ggsave("Fig.S11.tiff",scale = 1,
       plot = grid.arrange(plot.monkey.A, plot.monkey.B, ncol = 2),
       path = "C:/Users/weichunc/Desktop",
       width = 30, height = 35, units = "cm",dpi=320)



############################################# Human data ############################################################
# Human.obs:                                                                                                        #
# Measured pooled Plasma from 1977 - 2006 (Haugh et al., 2009)                                                      #          
# This is the only human study used in the calibration and optimization. Other datasets were used in the evaluation.#                                                                                                                   #
#####################################################################################################################
## Human plot

Sim.human  = as.data.frame(pred.human (Human.MCMC$bestpar[1:17], pred=TRUE))
human.sim  = cbind.data.frame (Time = Sim.human$Time, CA = Sim.human$CA)

plot.human =
  ggplot() +
  geom_line  (data = human.sim,  aes(x = Time, y = CA), color = "red", size=1.2) +
  geom_point (data = Human.obs,  aes(x = Time, y = CA), color = "black",size=2.5) 

# plot.human  =  plot.human  +
#   annotation_logticks(sides = "l") +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x)),
#                 limits = c(1e-1,1e2))+ xlim(c(0,35))

plot.human  = plot.human  + Theme.Fig + labs(x="",y="")
plot.human

### Save figure #####
ggsave("Fig.S12.tiff",scale = 1,
       plot = plot.human,
       path = "C:/Users/weichunc/Desktop",
       width = 35, height = 20, units = "cm",dpi=320)



