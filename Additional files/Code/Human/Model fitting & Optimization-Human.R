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
mod <- mcode ("pbpk", HumanPBK.code)


############################################# Model Calibration with levenberg-marquart #############
# Fitting and plot for humna study from norway                                                      #
# A1: Measured pooled Plasma from 1977 - 2006 (Haug et al., 2009)                                  #
# Concentration of PFCs (ng/ml) in pooled serum from men, age 40-50                                 #
#####################################################################################################

# input callibrated data set; (Haug et al., 2009) 
Human.obs <- read.csv(file="A1.csv")
names(Human.obs)=c("Time", "CA")

# Exposure scenario
BW          = 82.3
tinterval   = 24
TDoses      = 365*25


PDOSEoral1  = 0.003                       # ug/kg/day Oral dose; daily dose to 2.3e-3 - 3.0e-3 ug/kg/day in 1999-2000  (Loccisano et al., 2011)
DOSEoral1   = PDOSEoral1*BW
PDOSEoral2  = 0.0019                      # ug/kg/day Oral dose; daily dose to 1.6e-3 - 1.9e-3 ug/kg/day in 2000-2003 (Loccisano et al., 2011)
DOSEoral2   = PDOSEoral2*BW

# (PFOS Plasmas conc. for 1999-2000 with exposure to 2.3-3.0e-3 ug/kg/day)
# (PFOS Plasmas conc. for 2003-2004 with exposure to 1.6-1.9e-3 ug/kg/day)

oral1 <- ev(ID=1, time=0, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
oral2 <- ev(ID=1, time=24*365*25,amt= DOSEoral2, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
ex <- oral1+oral2 

# set up the exposure time
tsamp=tgrid(0,tinterval*(TDoses-1)+24*365*10,24*365)  # 25+10 = 35 years


# prediction function
Pred <- function(pars, pred=FALSE) {
  
  ##' Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  out <- 
    mod %>% 
    param(pars) %>%
    Req(Plasma)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid=tsamp)
    out<-cbind.data.frame (Time = out$time/(24*365),
                           CA   = out$Plasma)
  
}


## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out<-  Pred (pars)
  cost<- modCost(model=out, obs= Human.obs,x="Time")
  return(cost)
}


## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters

theta.int <- log(c(
  Vmax_baso_invitro          = 439.20,             ## Vmax of basolateral transpoters
  Km_baso                    = 20100,              ## Km of basolateral transpoters
  Vmax_apical_invitro        = 37400,              ## Vmax of apical transpoters
  Km_apical                  = 77500,              ## Km of apical transpoters
  RAFapi                     = 0.0007,             ## Relative acitivty factor of apical transpoter
  RAFbaso                    = 1,                  ## Relative activity factor of basolateral transpoter     
  KeffluxC                   = 0.100,              ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                     = 0.0001,             ## Biliary elimiatnion rate
  KurineC                    = 0.063,              ## Urinary elimination rate
  Free                       = 0.025,              ## Free fraction in plasma
  PL                         = 2.67,               ## Liver/plasma partition coefficient (PC)
  PK                         = 1.26,               ## Kidney/plasma PC
  PRest                      = 0.2,                ## Restofbody/plasma PC
  K0C                        = 1,                  ## Rate of absorption of PFOS in the stomach
  Kdif                       = 0.001,              ## Diffusion rate from proximal tubule cells;
  Kabsc                      = 2.120,              ## Rate of absorption of PFOS in the small intestines
  KunabsC                    = 7.06e-5             ## Rate of unabobded dose to appear in feces
))
# Senstivity function (FME)
# Check the senstive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen=summary(SnsPlasma)
write.csv(Sen,file="Sen.human.csv")

plot(summary(SnsPlasma))

## set up senstivie or necessary parametesr as model input
theta <- log(c(
  Vmax_baso_invitro          = 439.20,             ## Vmax of basolateral transpoters
  #Km_baso                   = 20100,              ## Km of basolateral transpoters
  Vmax_apical_invitro        = 37400,              ## Vmax of apical transpoters
  Km_apical                  = 77500,              ## Km of apical transpoters
  RAFapi                     = 0.0007,             ## Relative acitivty factor of apical transpoter
  #RAFbaso                   = 1,                  ## Relative activity factor of basolateral transpoter     
  KeffluxC                   = 0.100,              ## Rate of clearance of PFOS from proximal tubule cells into bloods
  KbileC                     = 0.0001,             ## Biliary elimiatnion rate
  KurineC                    = 0.063,              ## Urinary elimination rate
  Free                       = 0.025,              ## Free fraction in plasma
  PL                         = 2.67                ## Liver/plasma partition coefficient (PC)
  #PK                        = 1.26,               ## Kidney/plasma PC
  #PRest                     = 0.2,                ## Restofbody/plasma PC
  #K0C                       = 1,                  ## Rate of absorption of PFOS in the stomach
  #Kdif                      = 0.001,              ## Diffusion rate from proximal tubule cells;
  #Kabsc                     = 2.120,              ## Rate of absorption of PFOS in the small intestines
  #KunabsC                   = 7.06e-5             ## Rate of unabobded dose to appear in feces
))


# Model fitting
# fit.iv.A=nls.lm(par=theta.iv, fn = ModCal, pred=TRUE, control = nls.lm.control(nprint=1))
#lower <- log(as.numeric(exp(theta)/10)%>%setNames(names(theta)))
#upper <- log(as.numeric(exp(theta)*10)%>%setNames(names(theta)))

system.time(Fit<- modFit(f=MCcost, p=theta, method ="Nelder-Mead",
             control = nls.lm.control(nprint=1)))

summary(Fit)
exp(Fit$par)

res=MCcost(Fit$par)$residuals$res
sum(res^2)


## Plot using calibarated parameters 
df=cbind.data.frame (Time  = Pred (Fit$par)$Time, 
                     CA    = Pred (Fit$par)$CA)


ggplot() +
  geom_line(data  = df,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Human.obs,aes(Time,CA),size=2) + ylab("Concentration") 


############################################# Model Calibration with MCMC ###################################################
## One data sets was used in model evaluation                                                                               #
## A1: Measured pooled Plasma from 1977 - 2006 (Haugh et al., 2009)                                                         #
#############################################################################################################################

## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC <- log(c(
  Vmax_baso_invitro          = 4.791432e+02,       ## Vmax of basolateral transpoters/ Calibarated parameters
  Km_baso                    = 20100,              ## Km of basolateral transpoters
  Vmax_apical_invitro        = 5.180353e+04,       ## Vmax of apical transpoters/ Calibarated parameters
  Km_apical                  = 6.437290e+04 ,      ## Km of apical transpoters/ Calibarated parameters
  RAFapi                     = 1.000342e-03,       ## Relative acitivty factor of apical transpoter/ Calibarated parameters
  RAFbaso                    = 1,                  ## Relative activity factor of basolateral transpoter     
  KeffluxC                   = 1.496175e-01,       ## Rate of clearance of PFOS from proximal tubule cells into bloods/ Calibarated parameters
  KbileC                     = 1.348380e-04 ,      ## Biliary elimiatnion rate/ Calibarated parameters
  KurineC                    = 9.696675e-02,       ## Urinary elimination rate/ Calibarated parameters
  Free                       = 1.444482e-02,       ## Free fraction in plasma/ Calibarated parameters
  PL                         = 2.029106e+00,       ## Liver/plasma partition coefficient (PC)/ Calibarated parameters
  PK                         = 1.26,               ## Kidney/plasma PC
  PRest                      = 0.2,                ## Restofbody/plasma PC
  K0C                        = 1,                  ## Rate of absorption of PFOS in the stomach
  Kdif                       = 0.001,              ## Diffusion rate from proximal tubule cells;
  Kabsc                      = 2.120,              ## Rate of absorption of PFOS in the small intestines
  KunabsC                    = 7.06e-5,            ## Rate of unabobded dose to appear in feces
  sig2                       = 0.5,                ## Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  sig_Vmax.baso              = 0.5,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Km.baso                = 0.5,
  sig_Vmax.apical            = 0.5,
  sig_Km.apical              = 0.5,
  sig_RAFapi                 = 0.5, 
  sig_RAFbaso                = 0.5,
  sig_KeffluxC               = 0.5,
  sig_KbileC                 = 0.5,
  sig_KurineC                = 0.5,
  sig_Free                   = 0.3,
  sig_PL                     = 0.3,
  sig_PK                     = 0.3,
  sig_PRest                  = 0.3,
  sig_K0c                    = 0.3,
  sig_Kabsc                  = 0.3,
  sig_Kdif                   = 0.3,
  sig_KunabsC                = 0.3
))

  
which_sig <- grep("sig", names(theta.MCMC))

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars){
  
  ## Get out of log domain
  pars.data <- lapply(pars [-which_sig],exp)
  names(pars.data) <- names(pars.data)

  out <- 
    mod %>% 
    param(pars.data) %>%
    Req(Plasma)%>%
    update(atol = 1E-8,maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid=tsamp)
  
  outdf <- cbind.data.frame (Time=out$time/(24*365),
                             CA=out$Plasma)

  ## Data mutate with prediction
  outdf = outdf [which(outdf$Time %in% Human.obs$Time),]

  ## log.Predition 
  log.yhat      <-log(outdf$CA)

  ## log.Observed data
  log.y         <-log(Human.obs$CA)

  # The method of Maximum likelihood
  sig2            <- as.numeric((exp(pars[which_sig][1])))
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
  
  # Maximau likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)
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
    mod <- mcode ("pbpk", HumanPBK.code)
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
MC.H.1 = as.mcmc (MCMC[[1]]$pars)     # first  chain
MC.H.2 = as.mcmc (MCMC[[2]]$pars)     # second chain
MC.H.3 = as.mcmc (MCMC[[3]]$pars)     # third  chain
MC.H.4 = as.mcmc (MCMC[[4]]$pars)     # fourth chain
combinedchains = mcmc.list(MC.H.1,MC.H.2,MC.H.3,MC.H.4) ## combine all chains
gelman.diag (combinedchains)          # Gelman convergence diagnosis
heidel.diag (combinedchains)          # covergence diagnosis/Heidelberger and Welch's convergence diagnostic
gelman.plot (combinedchains)          # gelman plot

# Save the poseior parameters (95% CI)
quan.human = exp(summary(MC.H.1)$quantiles)  


# Trace plot using bayesplot
## Covergences plot
color_scheme_set("black")
mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC),
  size = 0.5,
  facet_args = list(nrow = 2)) +
  ggplot2::scale_color_brewer()


## output the MCMC results
write.csv(quan.Human,file="Human.summary_pos.csv")
write.csv(MC.H.1,file="Human.pos.csv")
saveRDS(MCMC[[1]],file ='Human.MCMC.rds')
saveRDS(combinedchains,file='Human.comb.rds')


## Plot using MCMC parameters

df=cbind.data.frame (Time  = Pred (MCMC.1$bestpar[-which_sig])$Time, 
                     CA    = Pred (MCMC.1$bestpar[-which_sig])$CA)


ggplot() +
  geom_line(data  = df,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Human.obs,aes(Time,CA),size=2) + ylab("Concentration") 



























