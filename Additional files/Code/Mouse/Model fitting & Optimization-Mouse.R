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
mod <- mcode ("micepbpk", micePBPK.code)


################################################################################################## 
# Mouse datasets for model calibration                                                           #
# Measured Plasma, kidney, liver, urine and feces as matrix (Chang et al., 2012)                 #
# A1: CD1 mice oral single dose of 1 mg/kg; Matrix: plasma                                       #
# A2: CD1 mice oral single dose of 1 mg/kg; Matrix: liver                                        #
# A3: CD1 mice oral single dose of 1 mg/kg; Matrix: kidney                                       #
# A4: CD1 mice oral single dose of 1 mg/kg; Matrix: urine                                        #
##################################################################################################

## The observed data adopoted from Chang et al., 2012; single oral dose of 1 mg/kg; 
## input data set

Obs.A1 <- read.csv(file ="A1.csv")   ## A1 dataset: oral single dose of 1 mg/kg; Matrix: plasma (Unit: ug/ml)
names(Obs.A1)=c("Time", "CA")

Obs.A2 <-read.csv(file  ="A2.csv")   ## A2 dataset: oral single dose of 1 mg/kg; Liver (Unit: ug/g)
names(Obs.A2)=c("Time", "CL")

Obs.A3 <- read.csv(file ="A3.csv")   ## A3 dataset: oral single dose of 1 mg/kg; Kidney (Unit: ug/g)
names(Obs.A3)=c("Time", "CK")

Obs.A4 <-read.csv(file  ="A4.csv")   ## A4 dataset: oral single dose of 1 mg/kg; Urine (Unit: % dose)
names(Obs.A4)=c("Time", "Curine")


## Define the prediction function (for model fitting using levenberg-marquart algorithm)
pred.mouse <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)

  ## Define the exposure scenario for oral exposue to 1 mg/kg-day
  
  BW           = 0.025                             ## kg, Mice body weight
  tinterval    = 24                                ## hr, Time interval
  TDoses       = 1                                 ## Dose times
  PDOSEoral    = 1                                 ## mg/kg-day, Single oral dose
  DOSEoral     = PDOSEoral*BW                      ## mg, amount of oral dose
  ex.oral<- ev(ID=1, amt= DOSEoral,                ## Set up the exposure regimes
            ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)

  ## Set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*200,1)     ## Simulation time 24*200 hours (200 days)
  
  ## Get a prediction
  ## Out: oral exposure to 1 mg/kg-d, matrix: Plasma, Liver, Kidney, urine and feces
  out.A <- 
    mod %>% 
    param(pars) %>%
    Req(Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-70,maxsteps = 50000) %>%
    mrgsim_d(data = ex.oral, tgrid=tsamp)
    out.A<-cbind.data.frame(Time=out.A$time/24, 
                            CA=out.A$Plasma, 
                            CL=out.A$Liver, 
                            CK= out.A$Kidney, 
                            Curine=(out.A$Aurine/DOSEoral)*100,
                            Cfeces=(out.A$Afeces/DOSEoral)*100) 
  return(out.A)
}



## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
theta.int <- log(c(
  Vmax_baso_invitro     = 393.45,                      ## Vmax of basolateral transpoters
  Km_baso               = 27.2,                        ## Km of basolateral transpoters
  Vmax_apical_invitro   = 9300,                        ## Vmax of apical transpoters
  Km_apical             = 52.3,                        ## Km of apical transpoters
  RAFapi                = 4.07,                        ## Relative acitivty factor of apical transpoter
  RAFbaso               = 3.99,                        ## Relative activity factor of basolateral transpoter
  KeffluxC              = 2.49,                        ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                = 0.004,                       ## Biliary elimiatnion rate
  KurineC               = 1.6,                         ## Urinary elimination rate
  Free                  = 0.09,                        ## Free fraction in plasma
  PL                    = 3.72,                        ## Liver/plasma partition coefficient (PC)
  PK                    = 0.8,                         ## Kidney/plasma PC
  PRest                 = 0.2,                         ## Restofbody/plasma PC
  K0C                   = 1,                           ## Rate of absorption of PFOS in the stomach
  Kabsc                 = 2.12,                        ## Rate of absorption of PFOS in the small intestines
  Kdif                  = 0.001,                       ## Diffusion rate from proximal tubule cells;
  KunabsC               = 7.05e-5                      ## Rate of unabobded dose to appear in feces
))

## Senstivity function (FME)
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = pred.mouse, parms = theta.int,varscale = 1)
Sen=summary(SnsPlasma)

plot(summary(SnsPlasma))


## Cost function (FME)
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out<- pred.mouse(pars)
  cost<- modCost(model=out,obs=Obs.A1,weight='std',x="Time")
  cost<- modCost(model=out,obs=Obs.A2,weight='std',x="Time",cost=cost)
  cost<- modCost(model=out,obs=Obs.A3,weight='std',x="Time",cost=cost)
  cost<- modCost(model=out,obs=Obs.A4,weight='std',x="Time",cost=cost)
  return(cost)
}


## set up senstivie or necessary parametesr as model input
## Fitting parmaers

theta <- log(c(
  #Vmax_baso_invitro     = 393.45,                      ## Vmax of basolateral transpoters
  #Km_baso               = 27.2,                        ## Km of basolateral transpoters
  Vmax_apical_invitro    = 9300,                        ## Vmax of apical transpoters
  #Km_apical             = 52.3,                        ## Km of apical transpoters
  RAFapi                 = 4.07,                        ## Relative acitivty factor of apical transpoter
  #RAFbaso               = 3.99,                        ## Relative activity factor of basolateral transpoter
  KeffluxC               = 2.49,                        ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                 = 0.004,                       ## Biliary elimiatnion rate
  #KurineC               = 1.6,                         ## Urinary elimination rate
  Free                   = 0.09,                        ## Free fraction in plasma
  #PL                    = 3.72,                        ## Liver/plasma partition coefficient (PC)
  #PK                    = 0.8,                         ## Kidney/plasma PC
  #PRest                 = 0.2,                         ## Restofbody/plasma PC
  #K0C                   = 1,                           ## Rate of absorption of PFOS in the stomach
  Kabsc                  = 2.12,                        ## Rate of absorption of PFOS in the small intestines
  Kdif                   = 0.001                        ## Diffusion rate from proximal tubule cells;
  #KunabsC               = 7.05e-5                      ## Rate of unabobded dose to appear in feces
))

## PBPK model fitting 
## least squres fit using levenberg-marquart (method "Marq") or Nelder-Mead algorithm

Fit<- modFit(f=MCcost, p=theta, method ="Marq",
             control = nls.lm.control(nprint=1))

summary(Fit)                           ## Summary of fit
exp(Fit$par)                           ## Get the arithmetic value out of the log domain

res=MCcost(Fit$par)$residuals$res      ## Check the residual for each time points
sum(res^2)                             ## Total residuals 


## Plot using calibarated parameters

Sim.fit.A = pred.mouse (par=Fit$par)
df.sim.A  = cbind.data.frame (Time=Sim.fit.A$Time, 
                              CA=Sim.fit.A$CA, 
                              CL=Sim.fit.A$CL,
                              CK=Sim.fit.A$CK,
                              Curine=Sim.fit.A$Curine)

plot.A1=
  ggplot() +
  geom_line(data  = df.sim.A,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1 ,aes(Time, CA),size=2.5) + ylab("Concentration") 


plot.A2=
  ggplot() +
  geom_line(data  = df.sim.A,aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2 ,aes(Time, CL),size=2.5) + ylab("Concentration") 


plot.A3=
  ggplot() +
  geom_line(data  = df.sim.A,aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3 ,aes(Time, CK),size=2.5) + ylab("Concentration") 


plot.A4=
  ggplot() +
  geom_line(data  = df.sim.A,aes(Time,Curine), col="firebrick", lwd=2)+
  geom_point(data = Obs.A4 ,aes(Time, Curine),size=2.5) + ylab("Concentration")


plot.A1
plot.A2
plot.A3
plot.A4

#######################################################################################################################                                                                                                                 #      
# Mouse datasets for model Evalution with MCMC  (Chang et al., 2012)                                                  #
# B1: CD1 mice oral single dose of 20 mg/kg; plasma                                                                   #
# B2: CD1 mice oral single dose of 20 mg/kg; liver                                                                    #
# B3: CD1 mice oral single dose of 20 mg/kg; kidney                                                                   #
# B4: CD1 mice oral single dose of 20 mg/kg; urine                                                                    #
#######################################################################################################################

## Datasets adopoted from Chang et al., 2012; single oral dose of 20 mg/kg; 
## Input data set 

Obs.B1 <- read.csv(file="B1.csv")  # B1 dataset: CD1 mice oral single dose of 20 mg/kg; matrix: plasma (Unit: ug/ml)
names(Obs.B1)=c("Time", "CA")

Obs.B2 <-read.csv(file="B2.csv")   # B2 dataset: CD1 mice oral single dose of 20 mg/kgPFOS; matrix: liver (Unit: ug/g)
names(Obs.B2)=c("Time", "CL")

Obs.B3 <- read.csv(file="B3.csv")  # B3 dataset: CD1 mice oral single dose of 20 mg/kgPFOS; matrix: kidney (Unit: ug/g)
names(Obs.B3)=c("Time", "CK")

Obs.B4 <-read.csv(file="B4.csv")   # B4 dataset: CD1 mice oral single dose of 20 mg/kgPFOS; matrix: urine (Unit: % dose)
names(Obs.B4)=c("Time", "Curine")


## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC<-log(c(
  Vmax_baso_invitro     = 393.45,                      ## Vmax of basolateral transpoters
  Km_baso               = 27.2,                        ## Km of basolateral transpoters
  Vmax_apical_invitro   = 4.185016e+03,                ## Vmax of apical transpoters/ calibrated parameters
  Km_apical             = 52.3,                        ## Km of apical transpoters
  RAFapi                = 2.812619e+00,                ## Relative acitivty factor of apical transpoter/ calibrated parameters
  RAFbaso               = 3.99,                        ## Relative activity factor of basolateral transpoter
  KeffluxC              = 5.603682e+00,                ## Rate of clearance of PFOS from proximal tubule cells into bloods/ calibrated parameters
  KbileC                = 3.935049e-04,                ## Biliary elimiatnion rate/ calibrated parameters
  KurineC               = 1.6,                         ## Urinary elimination rate
  Free                  = 1.971933e-02,                ## Free fraction in plasma/ calibrated parameters
  PL                    = 3.72,                        ## Liver/plasma partition coefficient (PC)
  PK                    = 0.8,                         ## Kidney/plasma PC
  PRest                 = 0.2,                         ## Restofbody/plasma PC
  K0C                   = 1,                           ## Rate of absorption of PFOS in the stomach
  Kabsc                 = 1.104189e+00,                ## Rate of absorption of PFOS in the small intestines/ calibrated parameters
  Kdif                  = 4.617048e-05,                ## Diffusion rate from proximal tubule cells/ calibrated parameters
  KunabsC               = 7.05e-5,                     ## Rate of unabobded dose to appear in feces  
  sig2                  = 0.5,                         ## Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of all parametesr is 0.3)
  sig_Vmax.baso         = 0.5,             ## default value of 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Km.baso           = 0.5,
  sig_Vmax.apical       = 0.5,
  sig_Km.apical         = 0.5,
  sig_RAFapi            = 0.5, 
  sig_RAFbaso           = 0.5,
  sig_KeffluxC          = 0.5,
  sig_KbileC            = 0.5,
  sig_KurineC           = 0.5,
  sig_Free              = 0.3,
  sig_PL                = 0.3,
  sig_PK                = 0.3,
  sig_PRet              = 0.3,
  sig_K0c               = 0.3,
  sig_Kabsc             = 0.3,
  sig_Kdif              = 0.3,
  sig_KunabsC           = 0.3
))

which_sig <- grep("sig", names(theta.MCMC))

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- lapply(pars [-which_sig],exp)

  
  # Exposure scenario for single oral dose of 20 mg/kg
  
  BW                    = 0.025                            ## kg, Mice body weight
  tinterval             = 24                               ## hr, Time interval
  TDoses                = 1                                ## Dose times
  PDOSEoral.B           = 20                               ## mg/kg-d, Single oral dose
  DOSEoral.B            = PDOSEoral.B*BW                   ## mg, amount of oral dose
  ex.oral.B<- ev(ID=1, 
              amt= DOSEoral.B, 
              ii=tinterval, 
              addl=TDoses-1, 
              cmt="AST", 
              replicate = FALSE)

  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*200,1)             ## simulation time 24*200 hr (200 days)
  
  ## Get a prediction
  out.B <- 
    mod %>%
    param(pars.data) %>%
    Req(Plasma, Liver, Kidney, Aurine, Afeces)%>%
    update(atol = 1E-8, maxsteps= 5000) %>%
    mrgsim_d(data = ex.oral.B, tgrid=tsamp)

   out.B = cbind.data.frame(Time   = out.B$time/24, 
                            CA     = out.B$Plasma, 
                            CL     = out.B$Liver,
                            CK     = out.B$Kidney, 
                            Curine = (out.B$Aurine/DOSEoral.B)*100)
  
  if (pred) return (out.B)
  
  out.B = outdf[which(out.B$Time %in% Obs.B1$Time),]
  
  ## log-transformed prediction
  log.yhat.CA     <- log(out.B$CA)
  log.yhat.CL     <- log(out.B$CL)
  log.yhat.CK     <- log(out.B$CK)
  log.yhat.Curine <- log(out.B$Curine)
  
  ## log-transformed experimental data
  log.y.CA        <- log(Obs.B1$CA)
  log.y.CL        <- log(Obs.B2$CL)
  log.y.CK        <- log(Obs.B3$CK)
  log.y.Curine    <- log(Obs.B4$Curine)
  
  ## The method of Maximum likelihood
  log.yhat        <- c(log.yhat.CA,log.yhat.CL,log.yhat.CK, log.yhat.Curine)
  log.y           <- c(log.y.CA,log.y.CL,log.y.CK, log.y.Curine)
  sig2            <- as.numeric((exp(pars[which_sig][1])))
  
  log_likelihood  <- -2*sum((dnorm (log.y,
                             mean = log.yhat,
                             sd   = sqrt(sig2), 
                             log=TRUE)))
                
  return(log_likelihood)
                
}

## Define the Prior distributions: either normal or lognormal distribution
## nomral distribution
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
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # prior distribution for population vraicne; sigma2
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
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
    mod <- mcode ("micepbpk", micePBPK.code)
    modMCMC(f             = mcmc.fun, 
            p             = theta.MCMC, 
            niter         = 500000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covrate matrix
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "weight" for the initial model variance
            ntrydr        = 2,                ## delayed Rejection
            burninlength  = 250000,           ## number of initial iterations to be removed from output.
            outputlength  = 50000)            ## number of output iterations           

    }
)

#end time
print(Sys.time()-strt)

stopCluster(cl)   

## Performance four chains to check the convergences
MC.mouse.1 = as.mcmc (MCMC[[1]]$pars) # first  chain
MC.mouse.2 = as.mcmc (MCMC[[2]]$pars) # second chain
MC.mouse.3 = as.mcmc (MCMC[[3]]$pars) # third  chain
MC.mouse.4 = as.mcmc (MCMC[[4]]$pars) # fourth chain

combinedchains = mcmc.list(MC.mouse.1,MC.mouse.2,MC.mouse.3,MC.mouse.4) ## combine all chains
gelman.diag (combinedchains)          # Gelman convergence diagnosis

## Save the poseior parameters (95% CI)
quan.mouse = exp(summary(MC.mouse.1)$quantiles)  

## Trace plot using bayesplot
## Covergences plot
color_scheme_set("black")
mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC[1:17]),
  size = 0.5,
  facet_args = list(nrow = 2)) +
  ggplot2::scale_color_brewer()



# output the MCMC results
write.csv(quan.mouse,file="mouse.summary_pos.csv")
write.csv(MC.mouse.1,file="mouse.pos.csv")
saveRDS(MCMC[[1]],file ='mouse.MCMC.rds')
saveRDS(combinedchains,file='mouse.comb.rds')


## Plot using MCMC parameters
## MOdel validation using positerir parametesr
Sim.fit.MCMC.B = mcmc.fun (par = MCMC[[1]]$bestpar,pred=TRUE)

df.sim.MCMC.B = cbind.data.frame (Time=Sim.fit.MCMC.B$Time, 
                                  CA=Sim.fit.MCMC.B$CA, 
                                  CL=Sim.fit.MCMC.B$CL,
                                  CK=Sim.fit.MCMC.B$CK,
                                  Curine=Sim.fit.MCMC.B$Curine)

plot.B1 =
  ggplot() +
  geom_line(data = df.sim.MCMC.B,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Obs.B1 ,aes(Time, CA),size=2.5) + ylab("Concentration") 


plot.B2 =
  ggplot() +
  geom_line(data = df.sim.MCMC.B,aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.B2 ,aes(Time, CL),size=2.5) + ylab("Concentration") 

plot.B3 =
  ggplot() +
  geom_line(data = df.sim.MCMC.B,aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.B3 ,aes(Time, CK),size=2.5) + ylab("Concentration") 

plot.B4=
  ggplot() +
  geom_line(data = df.sim.MCMC.B,aes(Time,Curine), col="firebrick", lwd=2)+
  geom_point(data = Obs.B4 ,aes(Time, Curine),size=2.5) + ylab("Concentration")

plot.B1
plot.B2
plot.B3
plot.B4






