## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(reshape2)   # melt function to reshape the table
library(tidyverse)  # Needed for the pipe %>% operator

## Loading human, rat, mouse, monkey MCMC data, from line 557 of the mrgModelfitting code
Human.MCMC  <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC    <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC  <- readRDS(file = "mouse.MCMC.rds")
Monkey.MCMC <- readRDS(file = "Monkey.MCMC.rds")

## loading the theta names
theta.names <- readRDS(file = "theta.names.rds")
which_sig   <- grep("sig", theta.names)

########################################## Table 3 ################################ 
# Median (2.5%, 97.5%) of posterior distribution for mouse, rat, monkey and human #
###################################################################################

## Prior parameters
## Population mean (u)
theta.Human  = readRDS (file="theta.Human.Rds")
theta.Monkey = readRDS (file="theta.Monkey.Rds")
theta.Mouse  = readRDS (file="theta.Mouse.Rds")
theta.Rat    = readRDS (file="theta.Rat.Rds")

## M value: Mean of population mean u
## S vaelu: SD of population mean u: set the CV equal to 50%
mean.mouse  = exp(theta.Mouse  [-which_sig]) # M-value for mouse
mean.rat    = exp(theta.Rat    [-which_sig]) # M-value for rat
mean.monkey = exp(theta.Monkey [-which_sig]) # M-value for monkey
mean.human  = exp(theta.Human  [-which_sig]) # M-value for human
mouse.sd    = mean.mouse*0.5                 # S-value for mouse; 
rat.sd      = mean.rat*0.5                   # S-value for Rat
monkey.sd   = mean.monkey*0.5                # S-value for monkey
human.sd    = mean.human*0.5                 # S-value for human

## 50%, 2.5% and 97.5% for prior population mean u
a.mouse  = qnorm(0.025, mean = mean.mouse, sd = mouse.sd)
b.mouse  = qnorm(0.975, mean = mean.mouse, sd = mouse.sd)
c.mouse  = qnorm(0.5,mean = mean.mouse, sd = mouse.sd)

a.rat    = qnorm(0.025, mean = mean.rat, sd = rat.sd)
b.rat    = qnorm(0.975, mean = mean.rat, sd = rat.sd)
c.rat  = qnorm(0.5,mean = mean.rat, sd = rat.sd)

a.monkey = qnorm(0.025, mean = mean.monkey, sd = monkey.sd)
b.monkey = qnorm(0.975, mean = mean.monkey, sd = monkey.sd)
c.monkey  = qnorm(0.5,mean = mean.monkey,  sd = monkey.sd)

a.human  = qnorm(0.025, mean = mean.human, sd = human.sd)
b.human  = qnorm(0.975, mean = mean.human, sd = human.sd)
c.human  = qnorm(0.5,mean = mean.human, sd = human.sd)

## Median (2.5%, 97.5%) for Posterior parameters of population mean
## Estimate the quantile of posterior parameters
quan.Human  = summary(as.mcmc(Human.MCMC$pars))$quantiles  
quan.Rat    = summary(as.mcmc(Rat.MCMC$pars))$quantiles    
quan.Mouse  = summary(as.mcmc(Mouse.MCMC$pars))$quantiles     
quan.Monkey = summary(as.mcmc(Monkey.MCMC$pars))$quantiles  

exp(quan.Human)
quan.Rat
quan.Mouse
Quan.Monkey
