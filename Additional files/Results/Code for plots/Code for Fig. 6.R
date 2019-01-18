########################################## Fig. 6 ################################ 
# circular plot of sensitivity analysis                                           #
###################################################################################
## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)      # Needed for the pipe %>% operator
library(mrgsolve)   # Needed to run the main PBPK code
library(reshape)    # melt function to reshape the table
library(ggplot2)    # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(grid)       # for plotting the figure
library(lattice)    # for plotting the figure

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

## Loading human, rat, mouse, monkey MCMC data
Human.MCMC        <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC          <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC        <- readRDS(file = "mouse.MCMC.rds")
Monkey.MCMC       <- readRDS(file = "Monkey.MCMC.rds")

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")
which_sig         <- grep("sig", theta.names)

## Sensitivity analysis
Pred <- function (pars.mouse,pars.rat,pars.monkey,pars.human){
  
  which_sig <- grep("sig", theta.names)
  ## Get out of log domain
  pars.mouse <- lapply(pars.mouse [-which_sig],exp) # pars.mouse is the entire parameter set; pars.mouse [-which_sig] means to keep parameters with "sig" only, and then do exp transformation, then reassign to pars.mouse
  pars.rat <- lapply(pars.rat [-which_sig],exp)
  pars.monkey <- lapply(pars.monkey [-which_sig],exp)
  pars.human <- lapply(pars.human [-which_sig],exp)
  
  ## Repeat dose exposure scenario: 
  
  BW.mouse          = 0.025                               ## mouse body weight
  BW.rat            = 0.3                                 ## rat body weight
  BW.monkey         = 3.5                                 ## monkey body weight
  BW.human          = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses.mouse      = 1                                   ## The number of dosing in the mouse
  TDoses.rat        = 98                                  ## The number of dosing in the rat
  TDoses.monkey     = 182                                 ## The number of dosing in the monkey
  TDoses.human      = 365*25                              ## The number of dosing in the human
  
  
  PDOSEoral.mouse   = 1                                   ## mg/kg; BW Oral dose, Change et al., 2012
  PDOSEoral.rat     = 1                                   ## mg/kg; BW Oral dose, Seacat et al., 2003
  PDOSEoral.monkey  = 0.75                                ## mg/kg; BW Oral dose, Seacat et al., 2002
  PDOSEoral.human   = 0.0045                              ## mg/kg; BW Oral dose, Olsen et al., 2003b
  
  DOSEoral.mouse    = PDOSEoral.mouse*BW.mouse            ## mg; amount of oral dose
  DOSEoral.rat      = PDOSEoral.rat*BW.rat                ## mg; amount of oral dose
  DOSEoral.monkey   = PDOSEoral.monkey*BW.monkey          ## mg; amount of oral dose
  DOSEoral.human    = PDOSEoral.human*BW.human            ## mg; amount of oral dose
  
  ex.mouse         <- ev(ID=1, amt= DOSEoral.mouse, ii=tinterval, addl=TDoses.mouse-1, cmt="AST", replicate = FALSE)
  ex.rat           <- ev(ID=1, amt= DOSEoral.rat, ii=tinterval, addl=TDoses.rat-1, cmt="AST", replicate = FALSE)
  ex.monkey        <- ev(ID=1, amt= DOSEoral.monkey, ii=tinterval, addl=TDoses.monkey-1, cmt="AST", replicate = FALSE)
  ex.human         <- ev(ID=1, amt= DOSEoral.human, ii=tinterval, addl=TDoses.human-1, cmt="AST", replicate = FALSE)
  
  
  ## set up the exposure time
  tsamp.mouse     = tgrid(0,tinterval*(TDoses.mouse-1)+24*365,24)          ## Simulated for 24*365 hours after dosing, but only obtained data at 24 h
  tsamp.rat       = tgrid(0,tinterval*(TDoses.rat-1)+24*100,24)            ## SD rat oral daily dose to 1 mg/kg for 98 days,
  tsamp.monkey    = tgrid(0,tinterval*(TDoses.monkey-1)+24*365,24)         ## 182 days and simulated for 24*365 hours after dosing
  tsamp.human     = tgrid(0,tinterval*(TDoses.human-1)+24*365*10,24*365)   ## Simulated for 25 + 10 years
  
  
  ## Get a prediction
  out.mouse <- 
    mod.mouse %>%
    param(pars.mouse) %>%
    Req(AUC_CA,AUC_CL,AUC_CK)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = ex.mouse, tgrid = tsamp.mouse )%>%
    filter(time!=0) # The code can produce time-dependent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  outdf.mouse = cbind.data.frame (Time       = out.mouse$time/24, 
                                  AUC_CA     = out.mouse$AUC_CA,
                                  AUC_CL     = out.mouse$AUC_CL,
                                  AUC_CK     = out.mouse$AUC_CK) 
  
  out.rat <- 
    mod.rat  %>%
    param(pars.rat) %>%
    Req(AUC_CA,AUC_CL,AUC_CK)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = ex.rat, tgrid = tsamp.rat)%>%
    filter(time!=0)
  
  outdf.rat = cbind.data.frame (Time       = out.rat$time/24, 
                                AUC_CA     = out.rat$AUC_CA,
                                AUC_CL     = out.rat$AUC_CL,
                                AUC_CK     = out.rat$AUC_CK) 
  
  out.monkey<- 
    mod.monkey  %>%
    param(pars.monkey) %>%
    Req(AUC_CA,AUC_CL,AUC_CK)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = ex.monkey, tgrid = tsamp.monkey)%>%
    filter(time!=0)
  
  outdf.monkey = cbind.data.frame (Time       = out.monkey$time/24, 
                                   AUC_CA     = out.monkey$AUC_CA,
                                   AUC_CL     = out.monkey$AUC_CL,
                                   AUC_CK     = out.monkey$AUC_CK) 
  
  out.human<- 
    mod.human  %>%
    param(pars.human) %>%
    Req(AUC_CA,AUC_CL,AUC_CK)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = ex.human, tgrid=tsamp.human)%>%
    filter(time!=0)
  
  outdf.human = cbind.data.frame (Time        = out.human$time/24, 
                                  AUC_CA     = out.human$AUC_CA,
                                  AUC_CL     = out.human$AUC_CL,
                                  AUC_CK     = out.human$AUC_CK) 
  
  
  
  return (list("outdf.mouse"  = outdf.mouse,
               "outdf.rat"    = outdf.rat,
               "outdf.monkey" = outdf.monkey,
               "outdf.human"  = outdf.human))
  
}

pars.mouse  = Mouse.MCMC$bestpar
pars.rat    = Rat.MCMC$bestpar
pars.monkey = Monkey.MCMC$bestpar
pars.human  = Human.MCMC$bestpar
R = Pred(pars.mouse,pars.rat,pars.monkey,pars.human)

## Create the matrix for normalized sensitivity coefficient data
NSC_mouse   = matrix(nrow=17,ncol=3)
NSC_rat     = matrix(nrow=17,ncol=3)
NSC_monkey  = matrix(nrow=17,ncol=3)
NSC_human   = matrix(nrow=17,ncol=3)


for (i in 1:17) {
  pars.mouse.new      <- log(c(exp(pars.mouse[i])*1.01,exp(pars.mouse[-i]))) # Each cycle, generate a new value of parameter i (e.g., 10.0a), and delete parameter i, so that you can proceed to the next parameter i+1
  pars.rat.new        <- log(c(exp(pars.rat[i])*1.01,exp(pars.rat[-i])))
  pars.monkey.new     <- log(c(exp(pars.monkey[i])*1.01,exp(pars.monkey[-i])))
  pars.human.new      <- log(c(exp(pars.human[i])*1.01,exp(pars.human[-i])))
  Rnew                <- Pred(pars.mouse.new,pars.rat.new,pars.monkey.new,pars.human.new)
  delta.P.mouse       <- exp(pars.mouse[i])/(exp(pars.mouse[i])*0.01) # Basically, the ratio is 100 for each parameter.
  delta.P.rat         <- exp(pars.rat[i])/(exp(pars.rat[i])*0.01)
  delta.P.monkey      <- exp(pars.monkey[i])/(exp(pars.monkey[i])*0.01)
  delta.P.human       <- exp(pars.human[i])/(exp(pars.human[i])*0.01)
  
  ## Estimated the AUC
  ## Mouse
  Mouse.AUC.CA.new   =  Rnew$outdf.mouse %>% filter (Time == 1) %>% select (AUC_CA) # AUC_CA, AUC_CL, and AUC_CK were defined in the main mrgsolve code
  Mouse.AUC.CA.ori   =  R$outdf.mouse %>% filter (Time == 1) %>% select (AUC_CA)
  Mouse.AUC.CL.new   =  Rnew$outdf.mouse %>% filter (Time == 1) %>% select (AUC_CL)
  Mouse.AUC.CL.ori   =  R$outdf.mouse%>% filter (Time == 1) %>% select (AUC_CL)
  Mouse.AUC.CK.new   =  Rnew$outdf.mouse%>% filter (Time == 1) %>% select (AUC_CK)
  Mouse.AUC.CK.ori   =  R$outdf.mouse%>% filter (Time == 1) %>% select (AUC_CK)
  
  delta.AUC.CA.mouse    =  Mouse.AUC.CA.new - Mouse.AUC.CA.ori
  delta.AUC.CL.mouse    =  Mouse.AUC.CL.new - Mouse.AUC.CL.ori
  delta.AUC.CK.mouse    =  Mouse.AUC.CK.new - Mouse.AUC.CK.ori
  
  ## Rat
  Rat.AUC.CA.new   =  Rnew$outdf.rat %>% filter (Time == 98) %>% select (AUC_CA)
  Rat.AUC.CA.ori   =  R$outdf.rat %>% filter (Time == 98) %>% select (AUC_CA)
  Rat.AUC.CL.new   =  Rnew$outdf.rat %>% filter (Time == 98) %>% select (AUC_CL)
  Rat.AUC.CL.ori   =  R$outdf.rat %>% filter (Time == 98) %>% select (AUC_CL)
  Rat.AUC.CK.new   =  Rnew$outdf.rat %>% filter (Time == 98) %>% select (AUC_CK)
  Rat.AUC.CK.ori   =  R$outdf.rat %>% filter (Time == 98) %>% select (AUC_CK)
  
  delta.AUC.CA.rat    =  Rat.AUC.CA.new - Rat.AUC.CA.ori
  delta.AUC.CL.rat    =  Rat.AUC.CL.new - Rat.AUC.CL.ori
  delta.AUC.CK.rat    =  Rat.AUC.CK.new - Rat.AUC.CK.ori
  
  ## Monkey
  Monkey.AUC.CA.new   =  Rnew$outdf.monkey %>% filter(Time == 182)%>% select (AUC_CA)
  Monkey.AUC.CA.ori   =  R$outdf.monkey %>% filter(Time == 182)%>% select (AUC_CA)
  Monkey.AUC.CL.new   =  Rnew$outdf.monkey %>% filter(Time == 182)%>% select (AUC_CL)
  Monkey.AUC.CL.ori   =  R$outdf.monkey %>% filter(Time == 182)%>% select (AUC_CL)
  Monkey.AUC.CK.new   =  Rnew$outdf.monkey %>% filter(Time == 182)%>% select (AUC_CK)
  Monkey.AUC.CK.ori   =  R$outdf.monkey %>% filter(Time == 182)%>% select (AUC_CK)
  
  delta.AUC.CA.monkey    =  Monkey.AUC.CA.new - Monkey.AUC.CA.ori
  delta.AUC.CL.monkey    =  Monkey.AUC.CL.new - Monkey.AUC.CL.ori
  delta.AUC.CK.monkey    =  Monkey.AUC.CK.new - Monkey.AUC.CK.ori
  
  ## Human
  Human.AUC.CA.new      =  Rnew$outdf.human %>% filter (Time == 25*365) %>% select (AUC_CA)
  Human.AUC.CA.ori      =  R$outdf.human %>% filter (Time == 25*365) %>% select (AUC_CA)
  Human.AUC.CL.new      =  Rnew$outdf.human %>% filter (Time == 25*365) %>% select (AUC_CL)
  Human.AUC.CL.ori      =  R$outdf.human %>% filter (Time == 25*365) %>% select (AUC_CL)
  Human.AUC.CK.new      =  Rnew$outdf.human %>% filter (Time == 25*365) %>% select (AUC_CK)
  Human.AUC.CK.ori      =  R$outdf.human %>% filter (Time == 25*365) %>% select (AUC_CK)
  
  delta.AUC.CA.human    =  Human.AUC.CA.new - Human.AUC.CA.ori
  delta.AUC.CL.human    =  Human.AUC.CL.new - Human.AUC.CL.ori
  delta.AUC.CK.human    =  Human.AUC.CK.new - Human.AUC.CK.ori
  
  
  NSC_mouse   [i, 1]   <- as.numeric((delta.AUC.CA.mouse/Mouse.AUC.CA.ori) * delta.P.mouse)
  NSC_mouse   [i, 2]   <- as.numeric((delta.AUC.CL.mouse/Mouse.AUC.CL.ori) * delta.P.mouse)
  NSC_mouse   [i, 3]   <- as.numeric((delta.AUC.CK.mouse/Mouse.AUC.CK.ori) * delta.P.mouse)
  NSC_rat     [i, 1]   <- as.numeric((delta.AUC.CA.rat/Rat.AUC.CA.ori) * delta.P.rat)
  NSC_rat     [i, 2]   <- as.numeric((delta.AUC.CL.rat/Rat.AUC.CL.ori) * delta.P.rat)
  NSC_rat     [i, 3]   <- as.numeric((delta.AUC.CK.rat/Rat.AUC.CK.ori) * delta.P.rat)
  NSC_monkey  [i, 1]   <- as.numeric((delta.AUC.CA.monkey/ Monkey.AUC.CA.ori) * delta.P.monkey)
  NSC_monkey  [i, 2]   <- as.numeric((delta.AUC.CL.monkey/ Monkey.AUC.CL.ori) * delta.P.monkey)
  NSC_monkey  [i, 3]   <- as.numeric((delta.AUC.CK.monkey/ Monkey.AUC.CK.ori) * delta.P.monkey)
  NSC_human   [i, 1]   <- as.numeric((delta.AUC.CA.human/Human.AUC.CA.ori) * delta.P.human)
  NSC_human   [i, 2]   <- as.numeric((delta.AUC.CL.human/Human.AUC.CL.ori) * delta.P.human)
  NSC_human   [i, 3]   <- as.numeric((delta.AUC.CK.human/Human.AUC.CK.ori) * delta.P.human)
  
}


colnames (NSC_mouse)  = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK") 
rownames(NSC_mouse)   = theta.names[1:17]
colnames (NSC_rat)    = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK") 
rownames(NSC_rat)     = theta.names[1:17]
colnames (NSC_monkey) = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK") 
rownames(NSC_monkey)  = theta.names[1:17]
colnames (NSC_human)  = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK") 
rownames(NSC_human)   = theta.names[1:17]


##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################

Circle.plot <- function (melt.data){ # melt.data is an argument of Circle.plot function.

# Set a number of 'empty bar' to add at the end of each group
  empty_bar=3
  to_add = data.frame(matrix(NA, empty_bar*nlevels(as.factor(melt.data$group)), ncol(melt.data)) )
  colnames(to_add) = colnames(melt.data)
  to_add$group=rep(levels(as.factor(melt.data$group)), each=empty_bar)
  melt.data=rbind(melt.data, to_add)
  melt.data=melt.data %>% arrange(group)
  melt.data$id=seq(1, nrow(melt.data)) # id is the number of rows. In total, there were 68 rows.

# Get the name and the y position of each label
 label_data=melt.data
 number_of_bar=nrow(label_data) # in total, there were 68 rows, 17 parameters * 4 species
 angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
 label_data$hjust<-ifelse( angle < -90, 1, 0)
 label_data$angle<-ifelse(angle < -90, angle+180, angle)

 # prepare a data frame for base lines
 base_data=melt.data %>%
 group_by(group) %>%
 summarize(start=min(id), end=max(id) - empty_bar) %>%
 rowwise() %>%
 mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
 grid_data = base_data
 grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
 grid_data$start = grid_data$start - 1
 grid_data=grid_data[-1,]

# Make the plot
 windowsFonts(Times=windowsFont("Times New Roman"))
 
 p.cir.plot <- 
 ggplot(melt.data, aes(x=as.factor(id), y = abs(value*100), fill = group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
 geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
  
# Add a val=80/60/40/20 lines. I do it at the beginning to make sure barplots are OVER it.
 geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
 geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
 geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
 geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
  
# Add text showing the value of each 80/60/40/20 lines
 annotate("text", x = rep(max(melt.data$id),4), y = c(20, 40, 60, 80), label = c("20%", "40%", "60%", "80%") , color="red", size=3 , angle=0, fontface="bold", hjust=1) +
 geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.9) +
 ylim(-100,120) +
 theme_minimal() +
 theme(
    legend.position         = "none",
    text                    = element_text (family = "Times"),
    panel.background        = element_blank (),
    plot.background         = element_blank (),
    axis.text               = element_blank(),
    axis.title              = element_blank(),
    panel.grid              = element_blank(),
    plot.margin             = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha = 1, size = 3.5, angle= label_data$angle, inherit.aes = FALSE) +
  
   
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=1.0 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size= 4, fontface="bold", inherit.aes = FALSE)
  
return (p.cir.plot)
}

###################### Fig. 8a; AUC of plasma ####################
melt.mouse.CA        = melt(NSC_mouse[,1]) # Initially, NSC_mouse is a dataset with two rows. melt can reshape it to a dataset with two columns.
melt.mouse.CA$group  = c("Mouse") # Add a third column of group with the values of Mouse for all rows.
melt.rat.CA          = melt(NSC_rat[,1])
melt.rat.CA$group    = c("Rat")
melt.monkey.CA       = melt(NSC_monkey[,1])
melt.monkey.CA$group = c("Monkey")
melt.human.CA        = melt(NSC_human[,1])
melt.human.CA$group  = c("Human")

melt.data.CA         = rbind (melt.mouse.CA,melt.rat.CA,melt.monkey.CA,melt.human.CA)
melt.data.CA$par     = rep(rownames(melt.mouse.CA),4) # Get the row names (parameter names), repeat four times, so that you can the parameter names for each of the four species

p1.CA                = Circle.plot (melt.data.CA%>%filter(abs(value)>0.01))
p1.CA
###################### Fig. 8b; AUC of liver ####################
melt.mouse.CL        = melt(NSC_mouse[,2])
melt.mouse.CL$group  = c("Mouse")
melt.rat.CL          = melt(NSC_rat[,2])
melt.rat.CL$group    = c("Rat")
melt.monkey.CL       = melt(NSC_monkey[,2])
melt.monkey.CL$group = c("Monkey")
melt.human.CL        = melt(NSC_human[,2])
melt.human.CL$group  = c("Human")

melt.data.CL         = rbind (melt.mouse.CL,melt.rat.CL,melt.monkey.CL,melt.human.CL)
melt.data.CL$par     = rep(rownames(melt.mouse.CA),4)
p1.CL                = Circle.plot (melt.data.CL%>%filter(abs(value)>0.01))
p1.CL

###################### Fig. 8c; AUC of Kidney ####################
melt.mouse.CK        = melt(NSC_mouse[,3])
melt.mouse.CK$group  = c("Mouse")
melt.rat.CK          = melt(NSC_rat[,3])
melt.rat.CK$group    = c("Rat")
melt.monkey.CK       = melt(NSC_monkey[,3])
melt.monkey.CK$group = c("Monkey")
melt.human.CK        = melt(NSC_human[,3])
melt.human.CK$group  = c("Human")

melt.data.CK         = rbind (melt.mouse.CK,melt.rat.CK,melt.monkey.CK,melt.human.CK)
melt.data.CK$par     = rep(rownames(melt.mouse.CA),4)

p1.CK                = Circle.plot (melt.data.CK%>%filter(abs(value)>0.01))

p1.CK 

grid.arrange(p1.CA,p1.CL,p1.CK, nrow = 1)

####### Save plot #######
ggsave("Fig.6.tiff",scale = 1,
       plot = grid.arrange(p1.CA,p1.CL,p1.CK, nrow = 1),
       path = "D:/Desktop",
       width = 45, height = 20, units = "cm",dpi=320)


dev.off()


