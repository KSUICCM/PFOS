########################################## Fig. 5 ################################ 
# Time-course profiles in mouse, rat, monkey and human                           #
##################################################################################

## Load libraries
library(mrgsolve) # Needed to run the main PBPK code
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2) # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr) # Needed for the pipe %>% operator
library(reshape) # melt function to reshape the table
library(scales) # for plotting the figure
library(gridExtra) # for plotting the figure
library(grid) # for plotting the figure
library(lattice) # for plotting the figure

## input PBPK model
micePBPK.code     <- readRDS (file = "micePBPK.RDS")
monkeyPBPK.code   <- readRDS (file = "monkeyPBPK.rds")
ratPBPK.code      <- readRDS (file = "ratPBPK.rds")
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

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)

## Population mean (u)
theta.Human       <- readRDS (file ="theta.Human.Rds")
theta.Monkey      <- readRDS (file ="theta.Monkey.Rds")
theta.Mouse       <- readRDS (file ="theta.Mouse.Rds")
theta.Rat         <- readRDS (file ="theta.Rat.Rds")

## Loading human, rat, mouse, monkey MCMC data
Human.MCMC        <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC          <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC        <- readRDS(file = "mouse.MCMC.rds")
Monkey.MCMC       <- readRDS(file = "Monkey.MCMC.rds")

## Theme
windowsFonts(Times=windowsFont("Times New Roman"))
Theme.Fig <-theme(
  legend.position         = "none",
  legend.title            = element_blank(),
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect(colour = "black", fill=NA, size=3),
  panel.background        = element_rect (fill="white"),
  panel.grid              = element_blank(),
  axis.text               = element_text (size   = 25, colour = "black", face = "bold"),    # tick labels along axes 
  axis.line               = element_line(size    = 0.8, colour = "black"),
  axis.title              = element_text (size   = 25, colour = "black", face = "bold"))   # label of axes
  
Theme.Fig2 <-theme(
  legend.position         = "right",
  legend.title            = element_blank(),
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect(colour = "black", fill=NA, size=3),
  panel.background        = element_rect (fill="white"),
  panel.grid              = element_blank(),
  axis.text               = element_text (size   = 25, colour = "black", face = "bold"),    # tick labels along axes 
  axis.line               = element_line(size    = 0.8, colour = "black"),
  axis.title              = element_text (size   = 25, colour = "black", face = "bold"))   # label of axes

############################################# Model Evalution ###########################                                                                                                                 #      
# Rat                                                                                   #                                                          
# - SD rat oral daily dose to 0, 0.03, 0.13, 0.34, 1.33 mg/kg-d for 28 days and 98 days #
# matrix: plasma and liver, data from Seacat et al. (2003)                              #                                                          
#########################################################################################

## Prediction function
pred.rat <- function (pars.rat) {

## Exposure scenario for rat oral daily exposue to 0.5, 2, 5, 10 ppm 
  
  BW.rat                 = 0.45                                      ## Rat body weight; from seacat et al. (2003)
  tinterval.rat          = 24                                        ## Time interval
  TDoses.rat             = 98                                        ## times, Dose times
  
  PDOSEoral.rat.a        = 0.03                                      ## mg/kg, BW Oral dose
  PDOSEoral.rat.b        = 0.13                                      ## mg/kg, BW Oral dose
  PDOSEoral.rat.c        = 0.34                                      ## mg/kg, BW Oral dose
  PDOSEoral.rat.d        = 1.33                                      ## mg/kg, BW Oral dose
  
  DOSEoral.rat.a         = PDOSEoral.rat.a*BW.rat                                    ## mg, amount of oral dose
  DOSEoral.rat.b         = PDOSEoral.rat.b*BW.rat                                    ## mg, amount of oral dose
  DOSEoral.rat.c         = PDOSEoral.rat.c*BW.rat                                    ## mg, amount of oral dose
  DOSEoral.rat.d         = PDOSEoral.rat.d*BW.rat                                    ## mg, amount of oral dose
  
  ex.rat.a               <- ev(ID=1, amt= DOSEoral.rat.a,
                                 ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  ex.rat.b               <- ev(ID=1, amt= DOSEoral.rat.b,
                                 ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  ex.rat.c               <- ev(ID=1, amt= DOSEoral.rat.c,
                                 ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  ex.rat.d               <- ev(ID=1, amt= DOSEoral.rat.d,
                                 ii=tinterval.rat, addl=TDoses.rat-1,cmt="AST",replicate = FALSE)
  
  tsamp.rat              = tgrid(0,tinterval.rat*(TDoses.rat-1)+24*200,1)   ## Siumuation 24*98 + 24*200 hours (298 days in total) Use the longer simulation time so that it captures the data at 28 days
  

  ## Rat ouput
  pars.rat %<>% lapply(exp)
  names(pars.rat) <- names(pars.rat)
  pars.rat        <- pars.rat[-which_sig ]
  
  outdf.rat.a <- 
    mod.rat %>% 
    param(pars.rat) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>% # 2000 is big enough to ensure stable simulations. rtol is the key factor that determines the speed. The smaller of rtol, the longer of the simulation time.
    mrgsim_d (data = ex.rat.a, tgrid = tsamp.rat)
  
  outdf.rat.b <- 
    mod.rat %>% 
    param(pars.rat) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.rat.b, tgrid = tsamp.rat)
  
  outdf.rat.c <- 
    mod.rat %>% 
    param(pars.rat) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.rat.c, tgrid = tsamp.rat)
  
  outdf.rat.d <- 
    mod.rat %>% 
    param(pars.rat) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.rat.d, tgrid = tsamp.rat)
  
  outdf.rat.a <- cbind.data.frame(Time   = outdf.rat.a$time/24, 
                                  CA     = outdf.rat.a$Plasma,                              
                                  CL     = outdf.rat.a$Liver)      
  
  
  outdf.rat.b <- cbind.data.frame(Time   = outdf.rat.b$time/24, 
                                  CA     = outdf.rat.b$Plasma,                              
                                  CL     = outdf.rat.b$Liver)      
  
  outdf.rat.c <- cbind.data.frame(Time   = outdf.rat.c$time/24, 
                                  CA     = outdf.rat.c$Plasma,                              
                                  CL     = outdf.rat.c$Liver)      
  
  outdf.rat.d <- cbind.data.frame(Time   = outdf.rat.d$time/24, 
                                  CA     = outdf.rat.d$Plasma,                              
                                  CL     = outdf.rat.d$Liver)      
  
  return (list("outdf.rat.a"   = outdf.rat.a, 
               "outdf.rat.b"   = outdf.rat.b,
               "outdf.rat.c"   = outdf.rat.c,
               "outdf.rat.d"   = outdf.rat.d))
  
}


Newtime.r   = pred.rat(theta.Rat)$outdf.rat.a$Time  # this is the new time variable, now it has been changed to sample per day.
nrwo.r = length (Newtime.r)

# Create the matrix 
MC.rat.a.CA    = matrix(nrow = nrwo.r, ncol = 5000) # nrow is the number of rows or the number of hours. 298 days, the first exposure is on day 0. So 297*24+ the title row and the 0 time point row. 
MC.rat.a.CL    = matrix(nrow = nrwo.r, ncol = 5000)
MC.rat.b.CA    = matrix(nrow = nrwo.r, ncol = 5000)
MC.rat.b.CL    = matrix(nrow = nrwo.r, ncol = 5000)
MC.rat.c.CA    = matrix(nrow = nrwo.r, ncol = 5000)
MC.rat.c.CL    = matrix(nrow = nrwo.r, ncol = 5000)
MC.rat.d.CA    = matrix(nrow = nrwo.r, ncol = 5000)
MC.rat.d.CL    = matrix(nrow = nrwo.r, ncol = 5000)

## Input paramters

for(i in 1:500){
  
  j = i *10  
  pars.rat             = Rat.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 5000 sets from 50000 total sets

  MCdata               = pred.rat (pars.rat)
  MC.rat.a.CA    [,i]  = MCdata $outdf.rat.a$CA    # this dataset contains 5000 columns of CA. outdf.rat.a is one of the output objects of the function.
  MC.rat.a.CL    [,i]  = MCdata $outdf.rat.a$CL
  MC.rat.b.CA    [,i]  = MCdata $outdf.rat.b$CA
  MC.rat.b.CL    [,i]  = MCdata $outdf.rat.b$CL
  MC.rat.c.CA    [,i]  = MCdata $outdf.rat.c$CA
  MC.rat.c.CL    [,i]  = MCdata $outdf.rat.c$CL
  MC.rat.d.CA    [,i]  = MCdata $outdf.rat.d$CA
  MC.rat.d.CL    [,i]  = MCdata $outdf.rat.d$CL
  
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}

M.rat.a.CA  = MC.rat.a.CA %>% apply(1,mean) # get the mean value of MC.rat.a.CA
SD.rat.a.CA = MC.rat.a.CA %>% apply(1,sd)
M.rat.b.CA  = MC.rat.b.CA %>% apply(1,mean)
SD.rat.b.CA = MC.rat.b.CA %>% apply(1,sd)
M.rat.c.CA  = MC.rat.c.CA %>% apply(1,mean)
SD.rat.c.CA = MC.rat.c.CA %>% apply(1,sd)
M.rat.d.CA  = MC.rat.d.CA %>% apply(1,mean)
SD.rat.d.CA = MC.rat.d.CA %>% apply(1,sd)
  
NewMC.rat.a.CA = cbind.data.frame(
                 Time = Newtime.r, Mean = M.rat.a.CA, SD = SD.rat.a.CA) # Create this object to plot the histogram. This dataset only contains predicted values. The observed values are below.

NewMC.rat.b.CA = cbind.data.frame(
                 Time = Newtime.r, Mean = M.rat.b.CA, SD = SD.rat.b.CA)

NewMC.rat.c.CA = cbind.data.frame(
                 Time = Newtime.r, Mean = M.rat.c.CA, SD = SD.rat.c.CA)

NewMC.rat.d.CA = cbind.data.frame(
                 Time = Newtime.r, Mean = M.rat.d.CA, SD = SD.rat.d.CA)


NewMC.rat.a.CA %<>% filter(Time == 98) # only obtain the predicted value at 98 days for the histogram.
NewMC.rat.a.CA$Dose = c("0.03")
NewMC.rat.b.CA %<>% filter(Time == 98)
NewMC.rat.b.CA$Dose = c("0.13")
NewMC.rat.c.CA %<>% filter(Time == 98)
NewMC.rat.c.CA$Dose = c("0.34")
NewMC.rat.d.CA %<>% filter(Time == 98)
NewMC.rat.d.CA$Dose = c("1.33")

Rat.CA <- rbind (NewMC.rat.a.CA, NewMC.rat.b.CA, NewMC.rat.c.CA, NewMC.rat.d.CA) # combine the data points.
Rat.CA$Matrix <- c("Pre.Plasma") # create a matrix column as the predicted plasma values.

# Liver
Newtime.r   = pred.rat (theta.Rat)$outdf.rat.a$Time
M.rat.a.CL  = MC.rat.a.CL %>% apply(1,mean)
SD.rat.a.CL = MC.rat.a.CL %>% apply(1,sd)
M.rat.b.CL  = MC.rat.b.CL %>% apply(1,mean)
SD.rat.b.CL = MC.rat.b.CL %>% apply(1,sd)
M.rat.c.CL  = MC.rat.c.CL %>% apply(1,mean)
SD.rat.c.CL = MC.rat.c.CL %>% apply(1,sd)
M.rat.d.CL  = MC.rat.d.CL %>% apply(1,mean)
SD.rat.d.CL = MC.rat.d.CL %>% apply(1,sd)

NewMC.rat.a.CL = cbind.data.frame(
  Time = Newtime.r, Mean = M.rat.a.CL, SD = SD.rat.a.CL)

NewMC.rat.b.CL = cbind.data.frame(
  Time = Newtime.r, Mean = M.rat.b.CL, SD = SD.rat.b.CL)

NewMC.rat.c.CL = cbind.data.frame(
  Time = Newtime.r, Mean = M.rat.c.CL, SD = SD.rat.c.CL)

NewMC.rat.d.CL = cbind.data.frame(
  Time = Newtime.r, Mean = M.rat.d.CL, SD = SD.rat.d.CL)


NewMC.rat.a.CL %<>% filter(Time == 98)
NewMC.rat.a.CL$Dose = c("0.03")
NewMC.rat.b.CL %<>% filter(Time == 98)
NewMC.rat.b.CL$Dose = c("0.13")
NewMC.rat.c.CL %<>% filter(Time == 98)
NewMC.rat.c.CL$Dose = c("0.34")
NewMC.rat.d.CL %<>% filter(Time == 98)
NewMC.rat.d.CL$Dose = c("1.33")

Rat.CL <- rbind (NewMC.rat.a.CL, NewMC.rat.b.CL, NewMC.rat.c.CL, NewMC.rat.d.CL)
Rat.CL$Matrix <- c("Pre.Liver")
Rat <- rbind (Rat.CA,Rat.CL)

## Observed data from Seacat et al., 2003
Obs.Time      = c(98)                          # days,Treatment 98 days
Obs.Dose      = c("0.03","0.13","0.34","1.33") # double quotation marks are not necessary, but it is OK to use it as it could be used as the x axis label.
Obs.rat.CA.M  = c(4.04,17.1,43.9,148)          # ug/ml, mean of PFOS concentration in Plasma for male rat after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CA.SD = c(0.8,1.22,4.9,14)             # ug/ml, sd of PFOS concentration in Plasma for male rat after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CL.M  = c(23.8,74,358,568)             # ug/g,  mean of PFOS conc. in liver for male rat after dosing 14 weeks
Obs.rat.CL.SD = c(3.5,6.2,29,107)
Obs.rat.CA    = cbind.data.frame (Time = Obs.Time, Mean = Obs.rat.CA.M, SD = Obs.rat.CA.SD, Dose = Obs.Dose, Matrix = c("Obs.Plasma")) # Create a matrix of obs.plasma just to categorize the data.
Obs.rat.CL    = cbind.data.frame (Time = Obs.Time, Mean = Obs.rat.CL.M, SD = Obs.rat.CL.SD, Dose = Obs.Dose,Matrix = c("Obs.Liver"))
Obs.rat       = rbind (Obs.rat.CA,Obs.rat.CL)
Comb.Rat = rbind (Rat, Obs.rat)
Comb.Rat$Dose = factor(Comb.Rat$Dose, levels = c("0.03","0.13","0.34","1.33"))
Comb.Rat$Matrix = factor(Comb.Rat$Matrix,levels=c("Pre.Plasma","Obs.Plasma","Pre.Liver","Obs.Liver"))

## Histogram plot
p1.r<- 
  ggplot(Comb.Rat, aes(x=as.factor(Dose), y = Mean, fill= as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat="identity", color="black", size = 1.2, 
           position=position_dodge()) + 
  scale_fill_manual(values=c("red","pink","#999999", "#E69F00"))+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p1.r = p1.r + 
       annotation_logticks(sides = "l",size=1.2) + 
       scale_y_log10(breaks=c(0,10,100,1000,10000),expand = c(0,0),limits = c(1,5000))

## Time course profiles plot
## Make the Observed data           
Obs.Time   = c(28,98)            # days,Treatment 28 and 98 days                                                                 
Obs.rat.CA.M  = c(41.8, 148)     # ug/ml, mean of PFOS concentration in Plasma for male rat after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CA.SD = c(7.9, 14)       # ug/ml, sd of PFOS concentration in Plasma for male monkey after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CL.M  = c(282, 568)      # ug/g, mean of PFOS concentration in liver for male monkey after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)
Obs.rat.CL.SD = c(45, 107)       # ug/g, sd of PFOS concentration in liver for male monkey after dosing to 14 weeks; Obtained from tale 2 (Seacat et al., 2003)  
Obs.rat.CA = cbind.data.frame (Time = Obs.Time, M = Obs.rat.CA.M, SD = Obs.rat.CA.SD)
Obs.rat.CL = cbind.data.frame (Time = Obs.Time, M = Obs.rat.CL.M, SD = Obs.rat.CL.SD)


MC.rat.CA.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.rat.d.CA, 1, function(y_est) c( # 1 indicates row; y_est indicates values from each row of MC.rat.d.CA; MC.rat.d.CA 7130 rows * 5000 columns
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.10, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.90, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.rat.CL.plot <- cbind(
   Time = Newtime.r, 
  as.data.frame(t(apply(MC.rat.d.CL, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)


p2.r.L <- 
  ggplot() + 
  geom_ribbon(data = MC.rat.CL.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="yellowgreen", alpha=0.3) +
  geom_ribbon(data = MC.rat.CL.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="green4", alpha = 0.3) +
  geom_line(data= MC.rat.CL.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "lightgreen") +
  geom_line(data= MC.rat.CL.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data=Obs.rat.CL, aes(x=Time, y= M), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Obs.rat.CL, aes(x=Time,ymin= M-SD, ymax = M+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05))

p2.r.p <- 
  ggplot() + 
  geom_ribbon(data = MC.rat.CA.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="red", alpha = 0.25) +
  geom_ribbon(data = MC.rat.CA.plot, aes(x = Time, ymin = ci_10, ymax = ci_90), 
              fill="red", alpha = 0.35) +
  geom_line(data= MC.rat.CA.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "black") +
  geom_line(data= MC.rat.CA.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data=Obs.rat.CA, aes(x=Time, y= M), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Obs.rat.CA, aes(x=Time,ymin= M-SD, ymax = M+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")


## Adjusted the x or y axis
p2.r.L <- p2.r.L +
  scale_y_continuous(limits = c(0,1100),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,220), expand = c(0,0))

p2.r.p <- p2.r.p +
  scale_y_continuous(limits = c(0,390),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,220),expand = c(0,0))

## Modified the theme
## plot
p2.r.L = p2.r.L + Theme.Fig + labs (x ="Time (days)", y="")
p2.r.p = p2.r.p + Theme.Fig + labs (x ="Time (days)", y="")
p1.r   = p1.r + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")

p2.r.L
p2.r.p


##### show plot #########
grid.arrange(p2.r.p, p2.r.L, p1.r,nrow = 1)

### Save figure #####
ggsave("Fig.5a.tiff",scale = 1,
       plot = grid.arrange(p2.r.p, p2.r.L, p1.r,nrow = 1),
       path = "D:/Desktop",
       width = 36, height = 12, units = "cm",dpi=320)


############################################# Model Evalution ###########################                                                                                                                 #      
# Cynomolgus Monkeys                                                                    #               #                                                          
# - oral daily dose to 0, 0.03, 0.15, 0.75 mg/kg-d for 182 days                         #
# matrix: plasma and liver, data from seacat et al. (2002)                              #                                                          
#########################################################################################
## Prediction function
pred.monkey <- function(pars.monkey) {
  
  ## Exposure scenario for rat oral daily exposue to 0.03, 0.15, 0.75 mg/kg-d 
  
  BW.monkey                 = 3.2
  tinterval.monkey          = 24                                        ## Time interval
  TDoses.monkey             = 182                                       ## Dose times 182 days
  
  PDOSEoral.monkey.a        = 0.03                                      ## mg/kg, BW Oral dose
  PDOSEoral.monkey.b        = 0.15                                      ## mg/kg, BW Oral dose
  PDOSEoral.monkey.c        = 0.75                                      ## mg/kg, BW Oral dose

  DOSEoral.monkey.a         = PDOSEoral.monkey.a*BW.monkey                                    ## mg, amount of oral dose
  DOSEoral.monkey.b         = PDOSEoral.monkey.b*BW.monkey                                    ## mg, amount of oral dose
  DOSEoral.monkey.c         = PDOSEoral.monkey.c*BW.monkey                                    ## mg, amount of oral dose

  ex.monkey.a               <- ev(ID=1, amt= DOSEoral.monkey.a,
                               ii=tinterval.monkey, addl=TDoses.monkey-1,cmt="AST",replicate = FALSE)
  ex.monkey.b               <- ev(ID=1, amt= DOSEoral.monkey.b,
                               ii=tinterval.monkey, addl=TDoses.monkey-1,cmt="AST",replicate = FALSE)
  ex.monkey.c               <- ev(ID=1, amt= DOSEoral.monkey.c,
                               ii=tinterval.monkey, addl=TDoses.monkey-1,cmt="AST",replicate = FALSE)

  tsamp.monkey              = tgrid(0,tinterval.monkey*(TDoses.monkey-1)+24*500,1)   ## Siumuation 24*182 + 24*500 hours (682 days)
  
  
  ## Monkey ouput
  pars.monkey %<>% lapply(exp)
  names(pars.monkey) <- names(pars.monkey)
  pars.monkey        <- pars.monkey[-which_sig ]
  
  outdf.monkey.a <- 
    mod.monkey %>% 
    param(pars.monkey) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.monkey.a, tgrid = tsamp.monkey)
  
  outdf.monkey.b <- 
    mod.monkey %>% 
    param(pars.monkey) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.monkey.b, tgrid = tsamp.monkey)
  
  outdf.monkey.c <- 
    mod.monkey %>% 
    param(pars.monkey) %>%
    Req (Plasma, Liver)%>%
    update(atol = 1E-5,rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim_d (data = ex.monkey.c, tgrid = tsamp.monkey)
  
 
  outdf.monkey.a <- cbind.data.frame(Time   = outdf.monkey.a$time/24, 
                                     CA     = outdf.monkey.a$Plasma,                              
                                     CL     = outdf.monkey.a$Liver)      
  
  
  outdf.monkey.b <- cbind.data.frame(Time   = outdf.monkey.b$time/24, 
                                     CA     = outdf.monkey.b$Plasma,                              
                                     CL     = outdf.monkey.b$Liver)      
  
  outdf.monkey.c <- cbind.data.frame(Time   = outdf.monkey.c$time/24, 
                                     CA     = outdf.monkey.c$Plasma,                              
                                     CL     = outdf.monkey.c$Liver)      
  
  
  return (list("outdf.monkey.a"   = outdf.monkey.a, 
               "outdf.monkey.b"   = outdf.monkey.b,
               "outdf.monkey.c"   = outdf.monkey.c))
  
}

Newtime.M      = pred.monkey(theta.Monkey)$outdf.monkey.a$Time
nrwo.M         = length (Newtime.M)

# Create the matrix 
MC.monkey.a.CA    = matrix(nrow = nrwo.M , ncol = 5000)
MC.monkey.a.CL    = matrix(nrow = nrwo.M , ncol = 5000)
MC.monkey.b.CA    = matrix(nrow = nrwo.M , ncol = 5000)
MC.monkey.b.CL    = matrix(nrow = nrwo.M , ncol = 5000)
MC.monkey.c.CA    = matrix(nrow = nrwo.M , ncol = 5000)
MC.monkey.c.CL    = matrix(nrow = nrwo.M , ncol = 5000)

## Input paramters

for(i in 1:5000){
  
  j = i *10  
  pars.monkey             = Monkey.MCMC$pars    [j,]
  
  MCdata               = pred.monkey (pars.monkey)
  MC.monkey.a.CA    [,i]  = MCdata $outdf.monkey.a$CA
  MC.monkey.a.CL    [,i]  = MCdata $outdf.monkey.a$CL
  MC.monkey.b.CA    [,i]  = MCdata $outdf.monkey.b$CA
  MC.monkey.b.CL    [,i]  = MCdata $outdf.monkey.b$CL
  MC.monkey.c.CA    [,i]  = MCdata $outdf.monkey.c$CA
  MC.monkey.c.CL    [,i]  = MCdata $outdf.monkey.c$CL

  cat("iteration = ", i , "\n")
}

# M.monkey.a.CA  = MC.monkey.a.CA %>% apply(1,mean)
# SD.monkey.a.CA = MC.monkey.a.CA %>% apply(1,sd)
# M.monkey.b.CA  = MC.monkey.b.CA %>% apply(1,mean)
# SD.monkey.b.CA = MC.monkey.b.CA %>% apply(1,sd)
# M.monkey.c.CA  = MC.monkey.c.CA %>% apply(1,mean)
# SD.monkey.c.CA = MC.monkey.c.CA %>% apply(1,sd)
# 
# NewMC.monkey.a.CA = cbind.data.frame(
#   Time = Newtime.M, Mean = M.monkey.a.CA, SD = SD.monkey.a.CA)
# 
# NewMC.monkey.b.CA = cbind.data.frame(
#   Time = Newtime.M, Mean = M.monkey.b.CA, SD = SD.monkey.b.CA)
# 
# NewMC.monkey.c.CA = cbind.data.frame(
#   Time = Newtime.M, Mean = M.monkey.c.CA, SD = SD.monkey.c.CA)
# 
# 
# NewMC.monkey.a.CA %<>% filter(Time == 182)
# NewMC.monkey.a.CA$Dose = c("0.03")
# NewMC.monkey.b.CA %<>% filter(Time == 182)
# NewMC.monkey.b.CA$Dose = c("0.15")
# NewMC.monkey.c.CA %<>% filter(Time == 182)
# NewMC.monkey.c.CA$Dose = c("0.75")
# 
# monkey.CA <- rbind (NewMC.monkey.a.CA, NewMC.monkey.b.CA, NewMC.monkey.c.CA)
# monkey.CA$Matrix <- c("Pre.Plasma")

# Liver
M.monkey.a.CL  = MC.monkey.a.CL %>% apply(1,mean)
SD.monkey.a.CL = MC.monkey.a.CL %>% apply(1,sd)
M.monkey.b.CL  = MC.monkey.b.CL %>% apply(1,mean)
SD.monkey.b.CL = MC.monkey.b.CL %>% apply(1,sd)
M.monkey.c.CL  = MC.monkey.c.CL %>% apply(1,mean)
SD.monkey.c.CL = MC.monkey.c.CL %>% apply(1,sd)

NewMC.monkey.a.CL = cbind.data.frame(
  Time = Newtime.M, Mean = M.monkey.a.CL, SD = SD.monkey.a.CL)

NewMC.monkey.b.CL = cbind.data.frame(
  Time = Newtime.M, Mean = M.monkey.b.CL, SD = SD.monkey.b.CL)

NewMC.monkey.c.CL = cbind.data.frame(
  Time = Newtime.M, Mean = M.monkey.c.CL, SD = SD.monkey.c.CL)


NewMC.monkey.a.CL %<>% filter(Time == 182)
NewMC.monkey.a.CL$Dose = c("0.03")
NewMC.monkey.b.CL %<>% filter(Time == 182)
NewMC.monkey.b.CL$Dose = c("0.15")
NewMC.monkey.c.CL %<>% filter(Time == 182)
NewMC.monkey.c.CL$Dose = c("0.75")

monkey.CL <- rbind (NewMC.monkey.a.CL, NewMC.monkey.b.CL, NewMC.monkey.c.CL)
monkey.CL$Matrix <- c("Pre.Liver")

monkey <- monkey.CL


Obs.Time.M        = c(182)                   # Treatment 182 days
Obs.Dose.M        = c("0.03","0.15","0.75")  # mg/kg-day, Dose group is "0.03", "0.15", "0.75" mg/kg-day (Seacat et al., 2002) 
#Obs.monkey.CA.M   = c(15.8,82.6,173)         # ug/ml, mean of PFOS concentration in Plasma for male monkey after dosing to 182 days; Obtained from tale 2 (Seacat et al., 2002)
#Obs.monkey.CA.SD  = c(1.4,25.2,37)           # ug/ml, sd of PFOS concentration in Plasma for male monkey after dosing to 182 days; Obtained from tale 2 (Seacat et al., 2002)
Obs.monkey.CL.M   = c(17.3,58.8,395)         # ug/g, mean of PFOS concentration in liver for male monkey after dosing to 182 days; Obtained from tale 2 (Seacat et al., 2002)  
Obs.monkey.CL.SD  = c(4.7,19.5,24)           # ug/g, sd of PFOS concentration in liver for male monkey after dosing to 182 days; Obtained from tale 2 (Seacat et al., 2002)  
#Obs.monkey.CA     = cbind.data.frame (Time = Obs.Time.M, Mean = Obs.monkey.CA.M , SD = Obs.monkey.CA.SD, Dose = Obs.Dose.M, Matrix = c("Obs.Plasma"))
Obs.monkey.CL     = cbind.data.frame (Time = Obs.Time.M, Mean = Obs.monkey.CL.M, SD = Obs.monkey.CL.SD, Dose = Obs.Dose.M,Matrix = c("Obs.Liver"))
Obs.monkey        = Obs.monkey.CL
Comb.monkey       = rbind (monkey, Obs.monkey)
Comb.monkey$Dose  = factor(Comb.monkey$Dose, levels = c("0.03","0.15","0.75"))
Comb.monkey$Matrix= factor(Comb.monkey$Matrix, levels = c("Pre.Liver","Obs.Liver"))

## Histogram PLOT
p1.M<- 
  ggplot(Comb.monkey, aes(x=as.factor(Dose), y = Mean, fill= as.factor(Matrix))) + 
  geom_bar(stat="identity", color="black", size = 1.2, 
           position=position_dodge()) + 
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))

p1.M = p1.M + 
       annotation_logticks(sides = "l",size=1.2) + 
       scale_y_log10(breaks=c(0,10,100,1000,10000),expand = c(0,0),limits = c(1,5000))

## Time-course plot

MC.monkey.CA.plot <- cbind(
  Time = Newtime.M, 
  as.data.frame(t(apply(MC.monkey.c.CA, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.10, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.90, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.monkey.CL.plot <- cbind(
  Time = Newtime.M , 
  as.data.frame(t(apply(MC.monkey.c.CL, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.10, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.90, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

Obs.monkey.CA<-Comb.monkey%>%filter(Dose == "0.75" & Matrix == "Obs.Plasma")
Obs.monkey.CL<-Comb.monkey%>%filter(Dose == "0.75" & Matrix == "Obs.Liver")

p2.monkey.L <- 
  ggplot() + 
  geom_ribbon(data = MC.monkey.CL.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="yellowgreen", alpha=0.3) +
  geom_ribbon(data = MC.monkey.CL.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="green4", alpha = 0.3) +
  geom_line(data= MC.monkey.CL.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "black") +
  geom_line(data= MC.monkey.CL.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data=Obs.monkey.CL, aes(x=Time, y= Mean), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Obs.monkey.CL, aes(x=Time,ymin= Mean-SD, ymax = Mean+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05))

p2.monkey.p <- 
  ggplot() + 
  geom_ribbon(data = MC.monkey.CA.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
               fill="red", alpha = 0.25) +
  geom_ribbon(data = MC.monkey.CA.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="red", alpha = 0.35) +
  geom_line(data= MC.monkey.CA.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "black") +
  geom_line(data= MC.monkey.CA.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data=Obs.monkey.CA, aes(x=Time, y= Mean), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Obs.monkey.CA, aes(x=Time,ymin= Mean-SD, ymax = Mean+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")


## Adjusted the x or y axis
p2.monkey.L <- p2.monkey.L +
  scale_y_continuous(limits = c(0,2250),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,650),expand = c(0,0))

# p2.monkey.p <- p2.monkey.p +
#   scale_y_continuous(limits = c(0,800),expand = c(0,0))+ 
#   scale_x_continuous(limits = c(0,650),expand = c(0,0))


## plot
p2.monkey.L = p2.monkey.L + Theme.Fig + labs (x ="Time (days)", y="")
p1.M        = p1.M + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")


##### show plot #####
grid.arrange(p2.monkey.L, p1.M, nrow = 1)

### Save figure #####
ggsave("Fig.5b.tiff",scale = 1,
       plot = grid.arrange(p2.monkey.L, p1.M, nrow = 1),
       path = "D:/Desktop",
       width = 24, height = 12, units = "cm",dpi=320)


############################################# Model Evalution ########################################################                                                                                                                 #      
# Human                                                                                                               #                                                          
# a- (PFOS Plasmas conc. for 2000-2001 with exposure to 2.2-4.5e-3 ug/kg/day); matrix: plasma; Olsen et al. (2003b)   #
# b- (PFOS Plasmas conc. for 2003-2004 with exposure to 8.5e-4-1.8e-3 ug/kg/day); matrix: Plasma; Olsen et al. (2008) #
# c- (Human tissue data with exposure to 1.8e-3 ug/kg/day (Fabrega et al., 2014); matrix: liver, kidney, plasma       #
# d- (PFOS liver conc. for 2000-2001 with exposure to 2.2-4.5e-3 ug/kg/day); matrix: liver; Olsen et al. (2003a)      #                                                  
######################################################################################################################


## Prediction function
pred.human <- function(pars.human) {
  
  ## Define the exposure scenarior for human
  BW.human                = 82.3
  tinterval.human         = 24
  TDoses.human.a          = 365*50                             # exposure scenario a:  Simulation time from 1950s (3M began to produce PFOS) to 2000 (PFOS was phased out)
  TDoses.human.b          = 365*30                             # exposure scenario b:  Simulaiton time from 2001 to 2011
  TDoses.human.c          = 365*50                             # exposure scenaior c&d:  simulation time for 50 years because the mean age of subject is 54 years
  
  PDOSEoral.human.a       = 0.0045                             # The estimated daily exposure to PFOS in 2000-2001 was 2.2 - 4.5e-3 ug/kg from Loccisano et al. (2011). The higher bound was used in the simulation because the exposure was high before 2000.
  DOSEoral.human.a        = PDOSEoral.human.a*BW.human
  PDOSEoral.human.b       = 0.0018                             # The estimated daily exposure to PFOS in 2006 was 8.5E-4 - 1.8e-3 ug/kg Loccisano et al. (2011). The lower bound was used in the simulation because the exposure was low after 2003.
  DOSEoral.human.b        = PDOSEoral.human.b*BW.human
  PDOSEoral.human.c       = 0.0018                             # Dietary exposure of PFOS was estimated to 1.8e-3 ug/kg body weight/day (Domingo et al.,2012)
  DOSEoral.human.c        = PDOSEoral.human.c*BW.human
  PDOSEoral.human.d       = 0.0045                             # The estimated daily exposure to PFOS in 2000-2001 was 2.2 - 4.5e-3 ug/kg from Loccisano et al. (2011). 
  DOSEoral.human.d        = PDOSEoral.human.d*BW.human
  
  oral1                   <- ev(ID=1, time = 0, amt= DOSEoral.human.a, ii=tinterval.human, addl = TDoses.human.a-1, cmt="AST", replicate = FALSE)
  oral2                   <- ev(ID=1, time = 24*365*50,amt= DOSEoral.human.b, ii=tinterval.human, addl = TDoses.human.b-1, cmt="AST", replicate = FALSE)
  ex.human.a              <- oral1 + oral2 # exposure scenario a is for datasets a and b, so b is incorporated into exposure scenario a.
  ex.human.c              <- ev(ID=1,amt= DOSEoral.human.c , ii=tinterval.human, addl = TDoses.human.c-1, cmt="AST", replicate = FALSE)
  ex.human.d              <- ev(ID=1,amt= DOSEoral.human.c , ii=tinterval.human, addl = TDoses.human.c-1, cmt="AST", replicate = FALSE)
  
  tsamp.human             = tgrid(0,tinterval.human*(TDoses.human.c -1)+24*365*30,24*365)  # here total simulation time is 50 + 30 years.
  
  
  ## Human ouput
  pars.human %<>% lapply(exp)
  pars.human  <- pars.human[-which_sig ]
  
  ## Scenario a: 2000  2001 with exposure to 4.5 e-3 ug/kg/day (Olsen et al., 2003a) # delete "?€?" if you see it.
  # exposure scenario a is for datasets a and b, so b is incorporated into exposure scenario a.
  outdf.human.a  <- 
    mod.human   %>% 
    param (pars.human) %>%
    Req (Plasma,Liver,Kidney) %>%
    update (atol= 1e-08, rtol= 1e-08, maxsteps = 5000) %>%
    mrgsim_d (data = ex.human.a, tgrid = tsamp.human)%>%
    filter(!duplicated(time)) # the new time variable is the original time divided by 24, by rounding to an integer value, some may be duplicated, so we have to exclude duplicated time points.
  
  outdf.human.a <- cbind.data.frame (Time = outdf.human.a$time/(24*365),
                                     CA   = outdf.human.a$Plasma,
                                     CL   = outdf.human.a$Liver,      
                                     CK   = outdf.human.a$Kidney)      
  
  
## Exposure senario c: Fabrega et al., (2014)
  outdf.human.c  <- 
    mod.human   %>% 
    param(pars.human) %>%
    Req(Plasma,Liver,Kidney) %>%
    update (atol= 1e-08, rtol= 1e-08, maxsteps = 5000) %>%
    mrgsim_d(data = ex.human.c, tgrid = tsamp.human)%>%
    filter(!duplicated(time))
  
  outdf.human.c <- cbind.data.frame (Time = outdf.human.c$time/(24*365),
                                    CA   = outdf.human.c$Plasma,
                                    CL   = outdf.human.c$Liver,      
                                    CK   = outdf.human.c$Kidney)      
  
  ## Exposure senario d: Human liver data (Olsen et al., 2003b)
  outdf.human.d  <- 
    mod.human   %>% 
    param (pars.human) %>%
    Req (Plasma,Liver) %>%
    update (atol= 1e-08, rtol= 1e-08, maxsteps = 5000) %>%
    mrgsim_d (data = ex.human.d, tgrid = tsamp.human) %>%
    filter(!duplicated(time))
  
    outdf.human.d <-cbind.data.frame (Time = outdf.human.d$time/(24*365),
                                      CA   = outdf.human.d$Plasma,
                                      CL   = outdf.human.d$Liver)
  
  
  
  return (list("outdf.a" = outdf.human.a, 
               "outdf.c" = outdf.human.c,
               "outdf.d" = outdf.human.d))
 }

Newtime.h       = pred.human(theta.Human)$outdf.a$Time
nrow.h          = length (Newtime.h)

## Create the matrix 
MC.human.a.CA    = matrix(nrow = nrow.h, ncol = 500) 
MC.human.a.CL    = matrix(nrow = nrow.h, ncol = 500)
MC.human.c.CA    = matrix(nrow = nrow.h, ncol = 500)
MC.human.c.CL    = matrix(nrow = nrow.h, ncol = 500)
MC.human.c.CK    = matrix(nrow = nrow.h, ncol = 500)
MC.human.d.CA    = matrix(nrow = nrow.h, ncol = 500)
MC.human.d.CL    = matrix(nrow = nrow.h, ncol = 500)


## Input paramters
for(i in 1:500){
  
  j = i *100  
  pars.human              = Human.MCMC$pars    [j,]
  
  MCdata                  = pred.human (pars.human)
  MC.human.a.CA     [,i]  = MCdata $outdf.a$CA
  MC.human.a.CL     [,i]  = MCdata $outdf.a$CL
  MC.human.c.CA     [,i]  = MCdata $outdf.c$CA
  MC.human.c.CL     [,i]  = MCdata $outdf.c$CL
  MC.human.c.CK     [,i]  = MCdata $outdf.c$CK
  MC.human.d.CA     [,i]  = MCdata $outdf.d$CA
  MC.human.d.CL     [,i]  = MCdata $outdf.d$CL
  
  cat("iteration = ", i , "\n")
}

M.human.a.CA    = MC.human.a.CA %>% apply(1,mean,na.rm=T) # if NA, then remove it. This is not necessary for this dataset.
SD.human.a.CA   = MC.human.a.CA %>% apply(1,sd,na.rm=T)
M.human.a.CL    = MC.human.a.CL %>% apply(1,mean,na.rm=T)
SD.human.a.CL   = MC.human.a.CL %>% apply(1,sd,na.rm=T)
M.human.c.CA    = MC.human.c.CA %>% apply(1,mean,na.rm=T)
SD.human.c.CA   = MC.human.c.CA %>% apply(1,sd,na.rm=T)
M.human.c.CL    = MC.human.c.CL %>% apply(1,mean,na.rm=T)
SD.human.c.CL   = MC.human.c.CL %>% apply(1,sd,na.rm=T)
M.human.c.CK    = MC.human.c.CK %>% apply(1,mean,na.rm=T)
SD.human.c.CK   = MC.human.c.CK %>% apply(1,sd,na.rm=T)
M.human.d.CA    = MC.human.d.CA %>% apply(1,mean,na.rm=T)
SD.human.d.CA   = MC.human.d.CA %>% apply(1,sd,na.rm=T)
M.human.d.CL    = MC.human.d.CL %>% apply(1,mean,na.rm=T)
SD.human.d.CL   = MC.human.d.CL %>% apply(1,sd,na.rm=T)

NewMC.human.a.CA = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.a.CA, SD = SD.human.a.CA)

NewMC.human.a.CL = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.a.CL, SD = SD.human.a.CL)

NewMC.human.c.CA = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.c.CA, SD = SD.human.c.CA)

NewMC.human.c.CL = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.c.CL, SD = SD.human.c.CL)

NewMC.human.c.CK = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.c.CK, SD = SD.human.c.CK)

NewMC.human.d.CA = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.d.CA, SD = SD.human.d.CA)

NewMC.human.d.CL = cbind.data.frame(
  Time = Newtime.h, Mean = M.human.d.CL, SD = SD.human.d.CL)


NewMC.human.a.CA %<>% filter(Time == 50) # 3M started to product PFOS in the early 1950s (~1951), sampling time in a was 2000-2001 (Olsen et al., 2003b) 
NewMC.human.a.CA$Dose = c("A")
NewMC.human.a.CL %<>% filter(Time == 50)
NewMC.human.a.CL$Dose = c("A")
NewMC.human.c.CA %<>% filter(Time == 50) # In the spanish study (F? brega et al., 2014), samples were collected from 54 years old subjects (median). After 5-8 decades, the concentrations have reached steady states, so concentrations at 50 years or 80 years are similar. 
NewMC.human.c.CA$Dose = c("C")
NewMC.human.c.CL %<>% filter(Time == 50)
NewMC.human.c.CL$Dose = c("C")
NewMC.human.c.CK %<>% filter(Time == 50)
NewMC.human.c.CK$Dose = c("C")
NewMC.human.d.CL %<>% filter(Time == 50)
NewMC.human.d.CL$Dose = c("D")

human.CA <- NewMC.human.c.CA
human.CA$Matrix <- c("Pre.Plasma")

human.CL <- NewMC.human.c.CL
human.CL$Matrix <- c("Pre.Liver")

human.CK <- NewMC.human.c.CK
human.CK$Matrix <- c("Pre.Kidney")

human <- rbind (human.CA,human.CL,human.CK)

# Data collected from table 2 in spanish study (F? brega et al., 2014)
Obs.Time.H       = c(50)     # In the spanish study (F? brega et al., 2014), samples were collected from a 28-86 years old subject (mean = 54 years old). After 5-8 decades, the concentrations have reached steady states, so concentrations at 50 years or 80 years are similar.   
Obs.Dose.H       = c("C")
Obs.human.CA.M   = c(13.6)   # ug/l, PFOS concentration in palsma    
Obs.human.CA.SD  = c(6.3)     
Obs.human.CL.M   = c(102)    # ng/g, PFOS concentration in liver
Obs.human.CL.SD  = c(123)
Obs.human.CK.M   = c(75.6)   # ng/g, PFOS concentration in kidney
Obs.human.CK.SD  = c(61.2)
Obs.human.CA     = cbind.data.frame (Time = Obs.Time.H, Mean = Obs.human.CA.M , SD = Obs.human.CA.SD, Dose = Obs.Dose.H, Matrix = c("Obs.Plasma"))
Obs.human.CL     = cbind.data.frame (Time = Obs.Time.H, Mean = Obs.human.CL.M,  SD = Obs.human.CL.SD, Dose = Obs.Dose.H, Matrix = c("Obs.Liver"))
Obs.human.CK     = cbind.data.frame (Time = Obs.Time.H, Mean = Obs.human.CK.M,  SD = Obs.human.CK.SD, Dose = Obs.Dose.H, Matrix = c("Obs.Kidney"))

Obs.human        = rbind (Obs.human.CA,Obs.human.CL,Obs.human.CK)
Comb.human       = rbind (human, Obs.human)
Comb.human$Matrix= factor(Comb.human$Matrix, levels = c("Pre.Plasma","Obs.Plasma","Pre.Liver","Obs.Liver","Pre.Kidney","Obs.Kidney"))
Comb.human$Tissue= factor(c("Plasma","Liver","Kidney","Plasma","Liver","Kidney"),levels = c("Plasma","Liver","Kidney"))

## Hitogram plot
p1.h <- 
     ggplot(Comb.human, aes(x=as.factor(Tissue), y = Mean, fill= as.factor(Matrix))) + 
     geom_bar(stat="identity", color="black", size = 1.2, position=position_dodge()) + 
     scale_fill_manual(values=c("red","pink","#999999", "#E69F00", "yellow","blue"))+
     geom_errorbar(aes(ymin=Mean, ymax=(Mean)+SD), width=.2, position=position_dodge(.9))


p1.h = p1.h + 
       annotation_logticks(sides = "l", size=1.2) + 
       scale_y_log10(breaks=c(0,10,100,1000,10000),expand = c(0,0),limits = c(1,5000))


### Time course plot for human
### Data collected from Figure 15 in Loccisano et al. (2011)
Obs.Time        = c(50,55)
Obs.human.CA.M  = c(36,15)
Obs.human.CA.SD = c(9, 6)
Obs.human.CA    = cbind.data.frame (Time = Obs.Time, M = Obs.human.CA.M, SD = Obs.human.CA.SD)

Obs.Time.d        = c(55)    # assumed exposre from early of 1950s (3M start to prdouce PFOS) to 2001 - 2006 
# Obs.human.CA.M.d  = c(24.38, 5.6)
# Obs.human.CA.SD.d = c(8.3, 1.9)
# Obs.human.CA.d    = cbind.data.frame (Time = Obs.Time.d, M = Obs.human.CA.M.d, SD = Obs.human.CA.SD.d)
Obs.human.CL.M.d  = c(19.2)  # ng/g; PFOS conc. in liver of male donors (avg. age is 50 years); data from Olsen et al., 2003a 
Obs.human.CL.SD.d = c(18)
Obs.human.CL.d    = cbind.data.frame (Time = Obs.Time.d, M = Obs.human.CL.M.d, SD = Obs.human.CL.SD.d)

MC.human.A.plot <- cbind(
  Time = Newtime.h, 
  as.data.frame(t(apply(MC.human.a.CA, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.10, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.90, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)



MC.human.D.plot  <- cbind(
  Time = Newtime.h , 
  as.data.frame(t(apply(MC.human.a.CL, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.10, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.90, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

p2.human.A <- 
  ggplot() + 
  geom_ribbon(data = MC.human.A.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="red", alpha=0.3) +
  geom_ribbon(data = MC.human.A.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="red", alpha = 0.3) +
  geom_line(data= MC.human.A.plot, aes(x = Time, y = median_est), 
            size = rel(1.2), colour = "black") +
  geom_line(data= MC.human.A.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Obs.human.CA, aes(x=Time, y= M), shape = 1, colour = "darkgreen", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Obs.human.CA, aes(x=Time,ymin= M-SD, ymax = M+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour = "black")
  

p2.human.D <- 
  ggplot() + 
  geom_ribbon (data = MC.human.D.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
               fill ="yellowgreen", alpha = 0.25) +
  geom_ribbon (data = MC.human.D.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
               fill ="green4", alpha = 0.35) +
  geom_line   (data   = MC.human.D.plot, aes(x = Time, y = median_est), 
               size   = rel(1.2), colour = "black") +
  geom_line   (data   = MC.human.D.plot, aes(x = Time, y = median_est), 
               colour = "black") +
  geom_point  (data  = Obs.human.CL.d, aes(x=Time, y= M), shape = 1, colour = "black", 
               fill  = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Obs.human.CL.d, aes(x=Time,ymin= M-SD, ymax = M+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="darkgreen")


## Adjusted the x or y axis
p2.human.A <- p2.human.A +
  scale_y_continuous(limits = c(0,85),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,65),expand = c(0,0))

p2.human.D <- p2.human.D +
  scale_y_continuous(limits = c(0,180),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,65),expand = c(0,0))

## Modified the theme
## plot
p2.human.A = p2.human.A + Theme.Fig + labs (x ="Times (years)", y="")
p2.human.D = p2.human.D + Theme.Fig + labs (x ="Times (years)", y="")
p1.h       = p1.h + Theme.Fig2 + labs (x ="Organs", y="")

p2.human.A
p2.human.D

#####
grid.arrange(p2.human.A, p2.human.D,p1.h, nrow = 1)



### Save figure ##
ggsave("Fig.5c.tiff",scale = 1,
       plot = grid.arrange(p2.human.A, p2.human.D,p1.h, nrow = 1),
       path = "D:/Desktop",
       width = 36, height = 12, units = "cm",dpi=320)


ggsave("Fig.5c2.tiff",scale = 1,
       plot = p1.h,
       path = "D:/Desktop",
       width = 10, height = 10, units = "cm",dpi=320)


dev.off()







