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
mod <- mcode ("Monkeypbpk", MonkeyPBPK.code)


########################################################################################## 
# Fitting and plot for singal iv dose in Cynomolgus Monkeys                              #
# Two datasets were measured by Chang et al., 2012                                       #
# A1: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: plasma                      #
# A2: Cynomolgus Monkeys/ single iv dose of 2 mg/kg, matrix: urine                       #
##########################################################################################

# input data set
Obs.A1 <- read.csv(file="A1.csv")   # A1 dataset: single iv dose of 2 mg/kg; matrix is plasma (Unit: ng/ml)
names(Obs.A1)=c("Time", "CA")

Obs.A2 <-read.csv(file="A2.csv")    # A2 dataset: single iv dose of 2 mg/kg; matrix is urine (unit: %dose)
names(Obs.A2)=c("Time", "Curine")


## Define the model calibration function to estimate the predict residuals 
## (for iv route; matrix: plasma, urine)

Pred.monkey <- function(pars) {
  
  ##' Get out of log domain
  pars %<>% lapply(exp)

  ## Define the exposure scenario for single iv dose of 2 mg/kg
  
  BW = 3.5                               ## kg, monkey body weight
  tinterval = 24                         ## hr, time interval
  TDoses = 1                             ## Dose times
  PDOSEiv.A = 2                          ## mg/kg, single iv dose 
  DOSEiv.A = PDOSEiv.A*BW                ## mg, amount of iv dose
  ex.iv.A <- ev(ID=1,                    ## set up exposure regimes
             amt= DOSEiv.A, 
             ii=tinterval, 
             addl=TDoses-1, 
             cmt="APlas_free", 
             replicate = FALSE)        
  

  ## set up the exposure time
  tsamp.iv = tgrid(0,tinterval*(TDoses-1)+24*200,1)   ## exposure time is 24*200 hours (200 days); 
  
  
  ##' Get a prediction
  out.A <- 
    mod %>% 
    param(pars) %>%
    Req(Plasma, Aurine)%>%
    update(atol = 1E-70,maxsteps=50000) %>%
    mrgsim_d(data = ex.iv.A, tgrid = tsamp.iv)
    out.A<-cbind.data.frame(Time = out.A$time/24,
                            CA   = out.A$Plasma, 
                            Curine = (out.A$Aurine/DOSEiv.A)*100)
  
    return (out.A)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
theta.int <- log(c(
  Vmax_baso_invitro      = 439.2,               ## Vmax of basolateral transpoters                  
  Km_baso                = 20.1,                ## Km of basolateral transpoters 
  Vmax_apical_invitro    = 37400,               ## Vmax of apical transpoters
  Km_apical              = 77.5,                ## Km of apical transpoters
  RAFapi                 = 0.0007,              ## Relative acitivty factor of apical transpoter
  RAFbaso                = 1,                   ## Relative activity factor of basolateral transpoter
  KeffluxC               = 0.1,                 ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                 = 0.0001,              ## Biliary elimiatnion rate
  KurineC                = 0.062,               ## Urinary elimination rate
  Free                   = 0.025,               ## Free fraction in plasma
  PL                     = 3.72,                ## Partion coefficient of liver tssiue (kim et al. 2006)
  PK                     = 0.8,                 ## Partion coefficient of kidney tssiue (kim et al. 2006)
  PRest                  = 0.2,                 ## Partion coefficient of residual tssiue (kim et al. 2006)
  K0C                    = 1,                   ## Rate of absorption of PFOS in the stomach
  Kdif                   = 0.001,               ## Diffusion rate from proximal tubule cells; 
  Kabsc                  = 2.12,                ## Rate of absorption of PFOS in the small intestines
  KunabsC                = 7.05e-5              ## Rate of unabobded dose to appear in feces
))


# Senstivity function (FME)
# Check the senstive parameters in the model
Sns <- sensFun(func = Pred.monkey, parms=theta.int,varscale = 1)
Sen=summary(Sns)
write.csv(Sen,file="Sen.monkey.csv")

plot(summary(Sns))

## Cost function (FME)
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out.A<- Pred.monkey(pars)
  cost<- modCost(model=out.A,obs=Obs.A1,x="Time")
  cost<- modCost(model=out.A,obs=Obs.A2,x="Time",cost=cost)
  return(cost)
}


## Set up senstivie or necessary parametesr as model input
## Fitting parmaers

theta <- log(c(
  #Vmax_baso_invitro      = 439.2,               ## Vmax of basolateral transpoters                  
  #Km_baso                = 20.1,                ## Km of basolateral transpoters 
  Vmax_apical_invitro     = 37400,               ## Vmax of apical transpoters
  Km_apical               = 77.5,                ## Km of apical transpoters
  RAFapi                  = 0.0007,              ## Relative acitivty factor of apical transpoter
  #RAFbaso                = 1,                   ## Relative activity factor of basolateral transpoter
  #KeffluxC               = 0.1,                 ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                  = 0.0001,              ## Biliary elimiatnion rate
  KurineC                 = 0.062,               ## Urinary elimination rate
  Free                    = 0.025,               ## Free fraction in plasma
  #PL                     = 3.72,                ## Partion coefficient of liver tssiue (kim et al. 2006)
  #PK                     = 0.8,                 ## Partion coefficient of kidney tssiue (kim et al. 2006)
  PRest                   = 0.2                  ## Partion coefficient of residual tssiue (kim et al. 2006)
  #K0C                    = 1,                   ## Rate of absorption of PFOS in the stomach
  #Kdif                   = 0.001                ## Diffusion rate from proximal tubule cells; 
  #Kabsc                  = 2.12                 ## Rate of absorption of PFOS in the small intestines
  #KunabsC                = 7.05e-5              ## Rate of unabobded dose to appear in feces
))

## PBPK model fitting 
## model fitting using levenberg-marquart (method "Marq") algorithm

Fit<- modFit(f=MCcost, p=theta, method ="Marq",
             control = nls.lm.control(nprint=1))

summary(Fit ) ## Summary of fit 
exp(Fit$par)  ## Get the arithmetic value out of the log domain

res=MCcost(Fit$par)$residuals$res ## Check the residual for each time points
sum(res^2)  ## Total residuals  


## plot using calibatrated parameters

Sim.fit.A = Pred.monkey (Fit$par)
df.sim.A  = cbind.data.frame (Time   = Sim.fit.A$Time, 
                              CA     = Sim.fit.A$CA,
                              Curine = Sim.fit.A$Curine)

plot.A1 =
  ggplot() +
  geom_line(data = df.sim.A,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1 ,aes(Time, CA),size=2.5) + ylab("Concentration") 


plot.A2 =
  ggplot() +
  geom_line(data = df.sim.A,aes(Time,Curine), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2 ,aes(Time, Curine)) + ylab("Concentration") + xlim(c(0,100))



plot.A1
plot.A2

########################################## Model Evalution with MCMC  ###########################################
## Oral.obs.A-C:  Cynomolgus Monkeys oral daily dose to 0.03,0.15,0.75 mg/kg for 182 days and monitored 1 year  #
## B1: 0.03 mg/kg-d                                                                                             #
## B2: 0.15 mg/kg-d                                                                                             #
## B3: 0.75 mg/kg-d                                                                                             #
## Matrix: plasma, data from Seacat et al. (2002)                                                               #                                                 
#################################################################################################################

## Input validation data set 

Obs.B1 <- read.csv (file="B1.csv")
names(Obs.B1) = c("Time", "CA")

Obs.B2 <- read.csv (file="B2.csv")
names(Obs.B2) = c("Time", "CA")

Obs.B3 <- read.csv (file="B3.csv")
names(Obs.B3) = c("Time", "CA")


## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC <- log(c(
  Vmax_baso_invitro      = 439.2,               ## Vmax of basolateral transpoters                  
  Km_baso                = 20.1,                ## Km of basolateral transpoters 
  Vmax_apical_invitro    = 7.697192e+04,        ## Vmax of apical transpoters/ calibarated parameters
  Km_apical              = 4.528324e+01,        ## Km of apical transpoters/ calibarated parameters
  RAFapi                 = 1.407897e-03,        ## Relative acitivty factor of apical transpoter/ calibarated parameters
  RAFbaso                = 1,                   ## Relative activity factor of basolateral transpoter
  KeffluxC               = 0.1,                 ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                 = 7.894670e-04,        ## Biliary elimiatnion rate/ calibarated parameters
  KurineC                = 9.289343e-02,        ## Urinary elimination rate/ calibarated parameters
  Free                   = 1.665475e-02,        ## Free fraction in plasma/ calibarated parameters
  PL                     = 3.72,                ## Partion coefficient of liver tssiue (kim et al. 2006)
  PK                     = 0.8,                 ## Partion coefficient of kidney tssiue (kim et al. 2006)
  PRest                  = 1.559218e-01,        ## Partion coefficient of residual tssiue (kim et al. 2006) / calibarated parameters
  K0C                    = 1,                   ## Rate of absorption of PFOS in the stomach
  Kdif                   = 0.001,               ## Diffusion rate from proximal tubule cells; 
  Kabsc                  = 2.12,                ## Rate of absorption of PFOS in the small intestines
  KunabsC                = 7.05e-5,             ## Rate of unabobded dose to appear in feces
  sig2                   = 0.5,                 ## Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of all parametesr is 0.3)
  sig_Vmax.baso           = 0.5,                ## default value of 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Km.baso             = 0.5,
  sig_Vmax.apical         = 0.5,
  sig_Km.apical           = 0.5,
  sig_RAFapi              = 0.5, 
  sig_RAFbaso             = 0.5,
  sig_KeffluxC            = 0.5,
  sig_KbileC              = 0.5,
  sig_KurineC             = 0.5,
  sig_Free                = 0.3,
  sig_PL                  = 0.3,
  sig_PK                  = 0.3,
  sig_PRest               = 0.3,
  sig_K0c                 = 0.3,
  sig_Kabsc               = 0.3,
  sig_Kdif                = 0.3,
  sig_KunabsC             = 0.3
))

which_sig <- grep("sig", names(theta.MCMC))

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- lapply(pars [-which_sig],exp)

  ## Repeat dose Exposure scenario: 
  ## monkeys were administered 0.03, 0.15, and 0.75 mg/kg PFOS orally for 183 days and monitored for a year after dosing
  
  BW.A        = 3.9                    ## Monkey body weight
  BW.B        = 3.3
  BW.C        = 3.2
  tinterval   = 24                     ## Time interval
  TDoses.rep  = 182                    ## Dose times 183 days
  
  PDOSEoral.A = 0.03                   ## mg/kg; BW Oral dose
  PDOSEoral.B = 0.15                   ## mg/kg; BW Oral dose
  PDOSEoral.C = 0.75                   ## mg/kg; BW Oral dose
  
  DOSEoral.rep.A = PDOSEoral.A*BW.A    ## mg; amount of oral dose
  DOSEoral.rep.B = PDOSEoral.B*BW.B    ## mg; amount of oral dose
  DOSEoral.rep.C = PDOSEoral.C*BW.C    ## mg; amount of oral dose
  
  ex.rep.A <- ev(ID=1, amt= DOSEoral.rep.A, ii=tinterval, addl=TDoses.rep-1, cmt="AST", replicate = FALSE)
  ex.rep.B <- ev(ID=1, amt= DOSEoral.rep.B, ii=tinterval, addl=TDoses.rep-1, cmt="AST", replicate = FALSE)
  ex.rep.C <- ev(ID=1, amt= DOSEoral.rep.C, ii=tinterval, addl=TDoses.rep-1, cmt="AST", replicate = FALSE)
  
  ## set up the exposure time
  tsamp.rep=tgrid(0,tinterval*(TDoses.rep-1)+24*565,24)  ## 183 days and simulated for 24*565 hours after dosing
  
  
  ## Get a prediction
  out.B1 <- 
    mod %>%
    param(pars.data) %>%
    Req(Plasma)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d(data = ex.rep.A, tgrid=tsamp.rep)%>%
    filter(time!=0)
  
  out.B1 = cbind.data.frame (Time = out.B1$time/24, 
                             CA   = out.B1$Plasma) 
                           
  out.B2 <- 
    mod %>%
    param(pars.data) %>%
    Req(Plasma)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d(data = ex.rep.B, tgrid=tsamp.rep)%>%
    filter(time!=0)
  
  out.B2 = cbind.data.frame (Time   = out.B2$time/24, 
                             CA     = out.B2$Plasma) 
  
  out.B3<- 
    mod %>%
    param(pars.data) %>%
    Req(Plasma)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d(data = ex.rep.C, tgrid=tsamp.rep)%>%
    filter(time!=0)
  
  out.B3 = cbind.data.frame (Time   = out.B3$time/24, 
                              CA    = out.B3$Plasma) 
  
  if (pred) return (list("out.B1" = out.B1, 
                         "out.B2" = out.B2,
                         "out.B3" = out.B3))
  
  out.B1 = out.B1[which(out.B1$Time %in% Obs.B1$Time),]
  out.B2 = out.B2[which(out.B2$Time %in% Obs.B2$Time),]
  out.B3 = out.B3[which(out.B3$Time %in% Obs.B3$Time),]
  
  ## log-transformed prediction
  log.yhat.B1      <- log(out.B1$CA)
  log.yhat.B2      <- log(out.B2$CA)
  log.yhat.B3      <- log(out.B3$CA)

  ## log-transformed experimental data
  log.y.B1         <- log(Obs.B1$CA)
  log.y.B2         <- log(Obs.B2$CA)
  log.y.B3         <- log(Obs.B3$CA)
  
  ## The method of Maximum likelihood: log likelihood
  log.yhat        <- c(log.yhat.B1,log.yhat.B2,log.yhat.B3)
  log.y           <- c(log.y.B1,log.y.B2,log.y.B3)
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
  sig  <- as.numeric (exp(pars[which_sig][2:18]))       # Coefficient of variation from population variance; sigmal0
  sig2 <- as.numeric (exp(pars[which_sig][1]))          # error variances from model esidual
  
  mean           = exp(theta.MCMC[-which_sig])
  CV             = 0.5                                  # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
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
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014   
  
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
    mod <- mcode ("Monkeypbpk", MonkeyPBPK.code)
    modMCMC(f             = mcmc.fun, 
            p             = theta.MCMC, 
            niter         = 500,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covrate matrix
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 250,           ## number of initial iterations to be removed from output.
            outputlength  = 50)            ## number of output iterations           

    }
)

#end time
print(Sys.time()-strt)

stopCluster(cl)   


## Performance four chains to check the convergences
MC.monkey.1 = as.mcmc (MCMC[[1]]$pars) # first  chain
MC.monkey.2 = as.mcmc (MCMC[[2]]$pars) # second chain
MC.monkey.3 = as.mcmc (MCMC[[3]]$pars) # third  chain
MC.monkey.4 = as.mcmc (MCMC[[4]]$pars) # fourth chain

combinedchains = mcmc.list(MC.monkey.1,MC.monkey.2,MC.monkey.3,MC.monkey.4) ## combine all chains
gelman.diag (combinedchains)        # Gelman convergence diagnosis

# Save the poseior parameters (95% CI)
quan.monkey = exp(summary(MC.monkey.1)$quantiles)  

# Trace plot using bayesplot
# Covergences plot
color_scheme_set("black")
mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC[1:17]),
  size = 0.5,
  facet_args = list(nrow = 6,ncol=3)) +
  ggplot2::scale_color_brewer()



# output the MCMC results
write.csv(quan.monkey,file="monkey.summary_pos.csv")
write.csv(MC.monkey.1,file="monkey.pos.csv")
saveRDS(MCMC[[1]],file ='monkey.MCMC.rds')
saveRDS(combinedchains,file='monkey.comb.rds')


## Plot using optimization parameters
Sim.mcmc.B1 = mcmc.fun(par = MCMC[[1]]$bestpar,pred=TRUE)
df.mcmc.B1  = cbind.data.frame (Time=Sim.mcmc.B1$out.B1$Time, CA=Sim.mcmc.B1$out.B1$CA)

Sim.mcmc.B2 = mcmc.fun(par = MCMC[[1]]$bestpar[-which_sig],pred=TRUE)
df.mcmc.B2  = cbind.data.frame (Time=Sim.mcmc.B2$out.B2$Time, CA=Sim.mcmc.B2$out.B2$CA)

Sim.mcmc.B3 = mcmc.fun(par = MCMC[[1]]$bestpar[-which_sig],pred=TRUE)
df.mcmc.B3  = cbind.data.frame (Time=Sim.mcmc.B3$out.B3$Time, CA=Sim.mcmc.B3$out.B3$CA)


plot.B1 = ggplot() +
          geom_line(data=df.mcmc.B1,aes(Time,CA), col="firebrick", lwd = 1.5)+
          geom_point(data=Obs.B1,aes(Time,CA)) + ylab("Concentration") 


plot.B2 = ggplot() +
          geom_line(data=df.mcmc.B2,aes(Time,CA), col="firebrick", lwd = 1.5)+
          geom_point(data=Obs.B2,aes(Time,CA)) + ylab("Concentration") 

plot.B3 = ggplot() +
          geom_line(data=df.mcmc.B3,aes(Time,CA), col="firebrick", lwd = 1.5)+
          geom_point(data=Obs.B3,aes(Time,CA)) + ylab("Concentration") 


plot.B1
plot.B2
plot.B3

