## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot

## Set working direction to the data files

## Build mrgsolve-based PBPK Model
mod <- mcode ("Ratpbpk", RatPBPK.code)


################################################################################################################################## 
# Model calibration and plot for 4 data sets                                                                                     #
# A1. Chang et al. (2012)   : SD rat oral single dose to 4.2 mg/kg, matrix: plasma                                        #                  
# A2. 3M unpublished data   : SD rat oral single dose to 2 mg/kg, matrix: plasma                                                 #
# A3. kim et al. (2016)     : SD rat oral single dose to 2 mg/kg, matrix: plasma                                       #
# A4. Johnson et al. (1979) : SD rat iv exposure to 4.2 mg/kg, matrix: urine                                                      #
##################################################################################################################################

## input data set for model calibration/ oral

Obs.A1 <-read.csv(file="A1.csv")        # A1 dataset: 4.2 mg/kg/single oral dose/matrix: plasma (Unit: ug/ml)/ Chang et al., (2012) 
names(Obs.A1) = c("Time", "CA")

Obs.A2 <-read.csv(file="A2.csv")        # A2 dataset: 2 mg/kg/single oral dose/matrix: plasma (Unit: ug/ml)/ Luccisano et al., (2012) from 3M unpublished data
names(Obs.A2) = c("Time", "CA")

Obs.A3 <-read.csv(file="A3.csv")        # A3 dataset: 2 mg/kg/single oral dose/matrix: plasma (Unit: ug/ml)/ Luccisano et al., (2012) from 3M unpublished data
names(Obs.A3) = c("Time", "CA")

## input data set for model calibration/ IV
Obs.A4 <- read.csv(file="A4.csv")       # A4 dataset: 4.2 mg/kg/single iv dose/matrix: urine (%Dose)/ Luccisano et al., (2012) from Johnson et al., 1979
names(Obs.A4) = c("Time", "Curine")

## Define the prediction function (for least squres fit using levenberg-marquart algorithm)
pred.rat <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain

  ## Define the three exposure scenario
  ## Exposure scenario for oral exposue to 4.2 mg/kg-d
  
  BW          = 0.3                                       # Rat body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEoral.A = 4.2                                       # Single oral dose from Luccisano et al. (2012)
  DOSEoral.A  = PDOSEoral.A*BW                            # Amount of oral dose
  
  # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
  ex.oral.A   <- ev (ID = 1, amt= DOSEoral.A, 
                    ii = tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)

  
  ## Exposure scenario for oral exposue to 2 mg/kg-d
  
  PDOSEoral.B = 2                                         # Single oral dose from Luccisano et al. (2012)
  DOSEoral.B  = PDOSEoral.B*BW                            # Amount of oral dose
  ex.oral.B   <- ev(ID = 1, amt= DOSEoral.B,              
                    ii = tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)

  ## Exposure scenario for IV exposue to 4.2 mg/kg-d
  
  PDOSEiv.A = 4.2                                         # Single IV dose from Luccisano et al. (2012)
  DOSEiv.A = PDOSEiv.A*BW                                 # Amount of IV dose
  IV.A <- ev(ID = 1, amt= DOSEiv.A, 
             ii = tinterval, addl=TDoses-1, cmt="APlas_free", replicate = FALSE)
  ex.iv.A <- IV.A
  
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*100,1)
  
  
  ## Get a prediction
  ## Out A: oral exposure to 4.2 mg/kg-d, matrix: Plasma, urine and feces
  out.A <- 
    mod %>%                                               # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(Plasma, Aurine, Afeces)%>%                        # select model output
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.A, tgrid=tsamp)               # Set up the simulation run
    out.A<-cbind.data.frame(Time=out.A$time/24, 
                             CA=out.A$Plasma, 
                             Curine=((out.A$Aurine)/DOSEoral.A)*100, 
                             Cfeces=((out.A$Afeces)/DOSEoral.A)*100)
   
  ## Out B: oral exposure to 2 mg/kg-d, matrix: Plasma, urine and feces
  out.B <- 
    mod %>% 
    param(pars) %>%
    Req(Plasma, Aurine,Afeces)%>%
    update(atol = 1E-70,maxsteps=50000) %>%
    mrgsim_d (data = ex.oral.B, tgrid=tsamp)
    out.B<-cbind.data.frame(Time=out.B$time/24, 
                              CA=out.B$Plasma, 
                              Curine=((out.B$Aurine)/DOSEoral.B)*100, 
                              Cfeces=((out.B$Afeces)/DOSEoral.B)*100)
  
  ## Outiv.A: IV exposure to 4.2 mg/kg-d, matrix: Plasma, urine 
  out.C <- 
    mod %>% 
    param(pars) %>%
    Req(Plasma,Aurine)%>%
    update(atol = 1E-70,maxsteps=50000) %>%
    mrgsim_d (data = ex.iv.A, tgrid=tsamp)
    out.C = cbind.data.frame(Time=out.C$time/24, 
                               CA=out.C$Plasma,
                               Curine=((out.C$Aurine)/DOSEiv.A)*100) 
    
  
    return(list("out.A"=out.A,
                "out.B"=out.B,
                "out.C"=out.C))
  
  }
    

## Cost fuction (FME) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out<-  pred.rat (pars)
  cost<- modCost  (model=out$out.A, obs= Obs.A1, x="Time")
  cost<- modCost  (model=out$out.B, obs= Obs.A2, x="Time",cost=cost)
  cost<- modCost  (model=out$out.B, obs= Obs.A3, x="Time",cost=cost)
  cost<- modCost  (model=out$out.C, obs= Obs.A4, x="Time",cost=cost)
  return(cost)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
theta.int <- log(c(
  Vmax_baso_invitro              = 393.45,                      ## Vmax of basolateral transpoters
  Km_baso                        = 27.2,                        ## Km of basolateral transpoters
  Vmax_apical_invitro            = 9300,                        ## Vmax of apical transpoters
  Km_apical                      = 52.3,                        ## Km of apical transpoters
  RAFapi                         = 4.07,                        ## Relative acitivty factor of apical transpoter
  RAFbaso                        = 3.99,                        ## Relative activity factor of basolateral transpoter
  KeffluxC                       = 2.49,                        ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                         = 0.004,                       ## Biliary elimiatnion rate
  KurineC                        = 1.6,                         ## Urinary elimination rate
  Free                           = 0.09,                        ## Free fraction in plasma
  PL                             = 3.72,                        ## Liver/plasma partition coefficient (PC)
  PK                             = 0.8,                         ## Kidney/plasma PC
  PRest                          = 0.2,                         ## Restofbody/plasma PC
  K0C                            = 1,                           ## Rate of absorption of PFOS in the stomach
  Kabsc                          = 2.12,                        ## Rate of absorption of PFOS in the small intestines
  Kdif                           = 0.001,                       ## Diffusion rate from proximal tubule cells;
  KunabsC                        = 7.05e-5                      ## Rate of unabobded dose to appear in feces
))

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)

Sen=summary(SnsPlasma)

plot(summary(SnsPlasma))

## set up senstivie or necessary parametesr as model input
## fitting parmaers

theta <- log(c(
  #Vmax_baso_invitro              = 393.45,                      ## Vmax of basolateral transpoters
  #Km_baso                        = 27.2,                        ## Km of basolateral transpoters
  Vmax_apical_invitro             = 9300,                        ## Vmax of apical transpoters
  Km_apical                       = 52.3,                        ## Km of apical transpoters
  RAFapi                          = 3.99,                        ## Relative acitivty factor of apical transpoter
  RAFbaso                         = 4.07,                        ## Relative activity factor of basolateral transpoter
  KeffluxC                        = 2.49,                        ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                          = 0.004,                       ## Biliary elimiatnion rate
  #KurineC                        = 1.6,                         ## Urinary elimination rate
  #Free                           = 0.09,                        ## Free fraction in plasma
  PL                              = 3.72,                        ## Liver/plasma partition coefficient (PC)
  #PK                             = 0.8,                         ## Kidney/plasma PC
  PRest                           = 0.2,                         ## Restofbody/plasma PC
  #K0C                            = 1,                           ## Rate of absorption of PFOS in the stomach
  #Kabsc                          = 2.12,                        ## Rate of absorption of PFOS in the small intestines
  Kdif                            = 0.001                        ## Diffusion rate from proximal tubule cells;
  #KunabsC                        = 7.05e-5                      ## Rate of unabobded dose to appear in feces
))


## PBPK model fitting 
## least squres fit using levenberg-marquart (method "Marq") algorithm

Fit<- modFit(f=MCcost, p=theta, method ="Marq",
             control = nls.lm.control(nprint=1))

summary(Fit)                                 ## Summary of fit 
exp(Fit$par)                                 ## Get the arithmetic value out of the log domain

res=MCcost(Fit$par)$residuals$res            ## Check the residual for each time points
sum(res^2)                                   ## Total residuals                       


## Model calibration plot using ggplot2 
Sim.fit.A = pred.rat (Fit$par)$out.A         ## Time-course PFOS concentration profiles using estimated parameters under exposure senaior A
Sim.fit.B = pred.rat (Fit$par)$out.B         ## Simulaiton of exposure scenaior B
Sim.fit.C = pred.rat (Fit$par)$out.C         ## Simulation of exposure scenario C

df.sim.A1 = cbind.data.frame (Time=Sim.fit.A$Time, CA=Sim.fit.A$CA)
                           

df.sim.A2 = cbind.data.frame (Time=Sim.fit.B$Time, CA=Sim.fit.B$CA)


df.sim.A3 = cbind.data.frame (Time=Sim.fit.B$Time, CA=Sim.fit.B$CA)


df.sim.A4 = cbind.data.frame (Time=Sim.fit.C$Time, Curine=Sim.fit.C$Curine)


## Plot
plot.A1=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1, aes(Time, CA),size=2.5) + ylab("Concentration") 


plot.A2=
  ggplot() +
  geom_line (data = df.sim.A2,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2, aes(Time, CA)) + ylab("Concentration") 


plot.A3=
  ggplot() +
  geom_line (data = df.sim.A3,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3 ,aes(Time, CA),size=2.5) + ylab("Concentration") 


plot.A4=
  ggplot() +
  geom_line (data = df.sim.A4,aes(Time,Curine), col="firebrick", lwd=2)+
  geom_point(data = Obs.A4 ,aes(Time, Curine),size=2.5) + ylab("Concentration") 


############################################# Optimization with MCMC Analaysis ####################################################
## Four data sets was used in model optimization with MCMC                                                                        #
## B1. 3M unpublished data  : SD rat oral single dose to 15 mg/kg,matrix: plasma                                                  #
## B2. Chang et al. (2012)  : SD rat oral single dose to 15 mg/kg,matrix: urine                                                   #
## B3. kim et al. (2016)    : SD rat iv single dose to 2 mg/kg,matrix: plasma                                                     #
## B4. 3M unpublished data  : SD rat oral daily dose to 1 mg/kg for 28 days; matrix: plasm                                        #
###################################################################################################################################

## lnput data sets
Obs.B1 <- read.csv(file ="B1.csv")      # B1 dataset from 3M unpublished data  
names(Obs.B1)=c("Time", "CA")

Obs.B2 <- read.csv(file ="B2.csv")      # B2 dataset from Chang et al. (2012)
names(Obs.B2)=c("Time", "Curine")

Obs.B3 <- read.csv(file ="B3.csv")      # B3 dataset from Kim et al., 2016    
names(Obs.B3)=c("Time", "CA")

Obs.B4 <- read.csv (file ="B4.csv")     # B4 datasets from 3M unpublished data 
names(Obs.B4)=c("Time", "CA")


## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)
theta.MCMC<-log(c(
  Vmax_baso_invitro              = 393.45,                      ## Vmax of basolateral transpoters
  Km_baso                        = 27.2,                        ## Km of basolateral transpoters
  Vmax_apical_invitro            = 1.808436e+03,                ## Vmax of apical transpoters
  Km_apical                      = 2.785610e+02,                ## Km of apical transpoters
  RAFapi                         = 1.903151e+00,                ## Relative acitivty factor of apical transpoter
  RAFbaso                        = 4.147229e+00,                ## Relative activity factor of basolateral transpoter
  KeffluxC                       = 2.098085e+00 ,               ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                         = 2.662026e-03,                ## Biliary elimiatnion rate
  KurineC                        = 1.6,                         ## Urinary elimination rate
  Free                           = 0.09,                        ## Free fraction in plasma
  PL                             = 3.664653e+00 ,               ## Liver/plasma partition coefficient (PC)
  PK                             = 0.8,                         ## Kidney/plasma PC
  PRest                          = 2.637778e-01,                ## Restofbody/plasma PC
  K0C                            = 1,                           ## Rate of absorption of PFOS in the stomach
  Kabsc                          = 2.12,                        ## Rate of absorption of PFOS in the small intestines
  Kdif                           = 5.147260e-04,                ## Diffusion rate from proximal tubule cells;
  KunabsC                        = 7.05e-5,                     ## Rate of unabobded dose to appear in feces
  sig2                           = 0.5,                         ## Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  sig_Vmax.baso                  = 0.5,               ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Km.baso                    = 0.5,
  sig_Vmax.apical                = 0.5,
  sig_Km.apical                  = 0.5,
  sig_RAFapi                     = 0.5, 
  sig_RAFbaso                    = 0.5,
  sig_KeffluxC                   = 0.5,
  sig_KbileC                     = 0.5,
  sig_KurineC                    = 0.5,
  sig_Free                       = 0.3,
  sig_PL                         = 0.3,
  sig_PK                         = 0.3,
  sig_PRest                      = 0.3,
  sig_K0c                        = 0.3,
  sig_Kabsc                      = 0.3,
  sig_Kdif                       = 0.3,
  sig_KunabsC                    = 0.3
))


which_sig <- grep("sig", names(theta.MCMC))           ## Get a list of parameters without "sig" character

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- pars [-which_sig] %>% lapply(exp)      ## Return a list of exp (parametrs) from log scale
  
  ## Exposure scenario for oral exposue to 15 mg/kg-d
  BW               = 0.3                              ## kg, Rat body weight                                      
  tinterval        = 24                               ## hr, Time interval                                 
  TDoses           = 1                                ## Dose times                                    
  PDOSEoral.C      = 15                               ## mg/kg-d, Single oral dose                           
  DOSEoral.C       = PDOSEoral.C*BW                   ## mg, amount of oral dose                   
  ex.oral.C        = ev(ID   = 1, 
                         amt  = DOSEoral.C, 
                         ii   = tinterval, 
                         addl = TDoses-1, 
                         cmt  = "AST", 
                         replicate = FALSE)

  ## Exposure scenario for IV exposue to 2 mg/kg-d
  
  PDOSEiv.B        = 2                                 ## mg/kg-d, Single iv dose                               
  DOSEiv.B         = PDOSEiv.B*BW                      ## mg, amount of iv dose              
  ex.iv.B          <- ev(ID =1, 
                      amt   = DOSEiv.B, 
                      ii    =tinterval, 
                      addl  =TDoses-1, 
                      cmt   ="APlas_free", 
                      replicate = FALSE)
  
  ## Repeat dose Exposure scenario (daily dose to 1 mg/kg for 28 days)
  
  TDoses.rep       = 28                                ## times, Dose times
  
  PDOSEoral        = 1                                 ## mg/kg, BW Oral dose
  DOSEoral.rep     = PDOSEoral*BW                      ## mg, amount of oral dose
  ex.rep           <- ev(ID =1, 
                      amt   = DOSEoral.rep, 
                      ii    = tinterval, 
                      addl  = TDoses.rep-1, 
                      cmt="AST", replicate = FALSE)
  
  
  # set up the exposure time
  tsamp = tgrid (0,tinterval*(TDoses-1)+24*100,1)       ## simulation 24*100 hr (100 days)

  ## Simulation of exposure scenaior A and B (oral single dose to 15 mg/kg)
  out.A <- 
    mod %>% 
    param (pars.data) %>%
    Req (Plasma, Aurine, Afeces) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%          ## atol: Absolute tolerance parameter
    mrgsim_d (data = ex.oral.C, tgrid=tsamp) %>%
    filter(time!=0)
  
  out.A = cbind.data.frame(Time=out.A$time/24, 
                             CA=out.A$Plasma, 
                             Curine=((out.A$Aurine)/DOSEoral.C)*100,
                             Cfeces=((out.A$Afeces)/DOSEoral.C)*100)
  
  ## Simulation of exposure scenaior C (IV exposue to 2 mg/kg-d)
  out.B <- 
    mod %>% 
    param (pars.data) %>%
    Req (Plasma,Aurine) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%          ## atol: Absolute tolerance parameter
    mrgsim_d (data = ex.iv.B, tgrid=tsamp) %>%
    filter(time!=0)
  
  out.B = cbind.data.frame(Time=out.B$time/24, 
                             CA=out.B$Plasma,
                             Curine=((out.B$Aurine)/DOSEiv.B)*100) 
  
  ## Simulation of exposure scenaior D (Daily exposue to 1 mg/kg-d for 28 days)
  out.C <-
    mod %>%
    param(pars.data) %>%
    Req(Plasma) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%          ## atol: Absolute tolerance parameter
    mrgsim_d(data=ex.rep, tgrid=tsamp) %>%
    filter(time!=0)

  out.C = cbind.data.frame(Time=out.C$time/24, CA=out.C$Plasma)

  if (pred) return(list("out.A" = out.A,     ## Exposure scenario A and B : oral dose/ 15 mg/kg
                        "out.B" = out.B,     ## Exposure scenario B       : iv dose/ 2 mg/kg
                        "out.C" = out.C))    ## Exposure scenario C       : repeated daily dose/1 mg/kg

  ## Get the same time-course profiles with experiment data from the simulation
  out.A = out.A[which(out.A$Time %in% Obs.B1$Time),]    
  out.B = out.B[which(out.B$Time %in% Obs.B3$Time),]
  out.C = out.C[which(out.C$Time %in% Obs.B4$Time),]
  
  ## log.Predition 
  log.yhat.oral.CA     <-log(out.A$CA)
  log.yhat.oral.Curine <-log(out.A$Curine)
  log.yhat.iv.CA       <-log(out.B$CA)
  log.yhat.rep.CA      <-log(out.C$CA)
  
  ## log.Observed data
  log.y.oral.CA         <-log(Obs.B1$CA)
  log.y.oral.Curine     <-log(Obs.B2$Curine)
  log.y.iv.CA           <-log(Obs.B3$CA)
  log.y.rep.CA          <-log(Obs.B4$CA)
  
  ## The method of Maximum likelihood
  log.yhat        <- c(log.yhat.oral.CA ,log.yhat.oral.Curine,log.yhat.iv.CA,log.yhat.rep.CA)
  log.y           <- c(log.y.oral.CA,log.y.oral.Curine ,log.y.iv.CA,log.y.rep.CA)
  
  sig2            <- as.numeric((exp(pars[which_sig][1]))) # Get the parameter of sig2 from the parameters
  
  log_likelihood  <- -2*sum ((dnorm (log.y,
                                     mean =log.yhat,
                                     sd = sqrt(sig2),
                                     log = TRUE)))
  return(log_likelihood)
  
}

## Define the Prior distributions: either normal or lognormal distribution
## nomral distribution
Prior <- function(pars) {
  
  ## Population level
  # The likelihood for population mean (parameters)
  pars.data = exp(pars[-which_sig])
  sig  <- as.numeric (exp(pars[which_sig][2:18]))                 # Coefficient of variation from population variance; sigmal0
  sig2 <- as.numeric (exp(pars[which_sig][1]))                    # error variances from model esidual
  
  mean           = exp(theta.MCMC[-which_sig])
  CV             = 0.5                                            # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  sd             = mean*CV
  
  # Calculate likelihoods of each parameters; P(u|M,S)
  prior_pars     = dtruncnorm(pars.data, 
                          a = qnorm(0.025, mean = mean, sd = sd), 
                          b = qnorm(0.975, mean = mean, sd = sd), 
                          mean = mean, sd = sd ) 
  
  # The likelihood for population variance; P(sigmal^2|sigmal0^2)
  CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
  CV.sig         = exp(theta.MCMC[which_sig])[2:18]               # Singmal0
  alpha          = (2+1)/(CU^2)                                   # Shape parametrer of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  beta           = (alpha-1)*CV.sig^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  
  # Calculate likelihoods of model error (sig2) and population variance (sig) parameters
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # Prior distribution for population vraicne; sigma2
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # Error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
  ## individual level; P(theta|u,sigmal^2)
  mean_i         = prior_pars
  sd_i           = sqrt(prior_sig)
  prior_pars_i   = dtruncnorm (prior_pars, 
                               a = qnorm(0.025, mean = mean_i, sd = sd_i), 
                               b = qnorm(0.975, mean = mean_i, sd = sd_i), 
                               mean = mean_i, sd = sd_i) 
  
  # log-transformed (log-likelihoods of each parameters)
  log.pri.pars   = log (prior_pars)
  log.pri.sig    = log (prior_sig)
  log.pri.pars.i = log (prior_pars_i)
  log.pri.sig2   = log (prior_sig2)
  
  # maximau likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)
  MLE =  -2*sum(log.pri.pars, log.pri.sig , log.pri.pars.i,log.pri.sig2)  
  
  return(MLE)
}

#################### MCMC simulation with parallel computing ############################

detectCores()                                ## check the cores
cl<- makeCluster(detectCores())              ## use all cores in our system     
registerDoParallel(cl)                       ## registers a cluster of all the cores on our system

# start time
strt<-Sys.time()

# parallel
system.time(
  MCMC <- foreach( i = 1:4, .packages = c('mrgsolve','magrittr','FME','truncnorm','EnvStats','invgamma','dplyr')) %dopar% {
    mod <- mcode ("Ratpbpk", RatPBPK.code)
    modMCMC(f             = mcmc.fun, 
            p             = theta.MCMC, 
            niter         = 500000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covariate matrix
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "Weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 250000,           ## number of initial iterations to be removed from output.
            outputlength  = 50000)            ## number of output iterations           
  }
)

#end time
print(Sys.time()-strt)

stopCluster(cl)   

## Performance four chains to check the convergences
MC.rat.1 = as.mcmc (MCMC[[1]]$pars)   # first  chain
MC.rat.2 = as.mcmc (MCMC[[2]]$pars)   # second chain
MC.rat.3 = as.mcmc (MCMC[[3]]$pars)   # third  chain
MC.rat.4 = as.mcmc (MCMC[[4]]$pars)   # fourth chain
combinedchains = mcmc.list(MC.rat.1,MC.rat.2,MC.rat.3,MC.rat.4) ## combine all chains
gelman.diag (combinedchains)          # Gelman convergence diagnosis

## Save the poseior parameters (95% CI)
quan.rat = exp(summary(MC.rat.1)$quantiles)  

## Trace plot using bayesplot
## Covergences plot
mcmc_trace (
  combinedchains,
  pars       = names(theta.MCMC[1:17]),
  size       = 0.5,
  facet_args = list(nrow = 2)) +
  ggplot2::scale_color_brewer()



## output the MCMC results with .csv and .rds
write.csv(quan.rat,file="rat.summary_pos.csv")
write.csv(MC.rat.1,file="rat.pos.csv")
saveRDS(MCMC[[1]],file ='rat.MCMC.rds')
saveRDS(combinedchains,file='rat.comb.rds')


## Make the plot using MCMC parameters
## single oral and iv
Sim.mcmc.A    = mcmc.fun  (par=MCMC[[1]]$bestpar,pred=TRUE)$out.A
Sim.mcmc.B    = mcmc.fun  (par=MCMC[[1]]$bestpar,pred=TRUE)$out.B
Sim.mcmc.C    = mcmc.fun  (par=MCMC[[1]]$bestpar,pred=TRUE)$out.C

df.sim.mcmc.B1    = cbind.data.frame (Time=Sim.mcmc.A$Time,CA=Sim.mcmc.A$CA)

df.sim.mcmc.B2    = cbind.data.frame (Time=Sim.mcmc.A$Time,Curine=Sim.mcmc.A$Curine)

df.sim.mcmc.B3    = cbind.data.frame (Time=Sim.mcmc.B$Time,CA=Sim.mcmc.B$CA)

df.sim.mcmc.B4    = cbind.data.frame (Time=Sim.mcmc.C$Time,CA=Sim.mcmc.C$CA)

Plot.B1 =
  ggplot() +
  geom_line(data=df.sim.mcmc.B1,aes(Time,CA), col="darkblue", lwd=1)+
  geom_point(data=Obs.B1,aes(Time,CA)) +  
  ylab("Concentration")


Plot.B2 =
  ggplot() +
  geom_line(data=df.sim.mcmc.B2,aes(Time,Curine), col="darkblue", lwd=1)+
  geom_point(data=Obs.B2,aes(Time,Curine)) +  
  ylab("Concentration")


Plot.B3 =
  ggplot() +
  geom_line(data=df.sim.mcmc.B3,aes(Time,CA), col="darkblue", lwd=1)+
  geom_point(data=Obs.B3,aes(Time,CA)) +  
  ylab("Concentration")


Plot.B4 =
  ggplot() +
  geom_line(data=df.sim.mcmc.B4,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data=Obs.B4,aes(Time,CA)) + ylab("Concentration") 


Plot.B1
Plot.B2
Plot.B3
Plot.B4







