
########################################## Fig. 3 ################################ 
# Densities plot of prior and posterior parameters uncertainty distribution       #
###################################################################################
## loading R packages
library(magrittr)  # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)     # Needed for the pipe %>% operator
library(mrgsolve)  # Needed to run the main PBPK code
library(reshape)   # melt function to reshape the table
library(tidyverse) # Needed for the pipe %>% operator
library(ggjoy)     # Used to create Figure 3
library(ggplot2)   # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.


## Loading human, rat, mouse, monkey MCMC data
Human.MCMC        <- readRDS(file = "Human.MCMC.rds")
Rat.MCMC          <- readRDS(file = "Rat.MCMC.rds")
Mouse.MCMC        <- readRDS(file = "Mouse.MCMC.rds")
Monkey.MCMC       <- readRDS(file = "Monkey.MCMC.rds")

## loading the theta names
theta.names       <- readRDS(file = "theta.names.rds")

## Sampling from posterior parameters to generate the posterior distributions 
## Human posteior distributions
M.Human  <- exp(Human.MCMC$pars) %>% apply(2,mean) # Get the mean of each column. 1 indicates row; 2 indicates column. Each column has 50000 values.

SD.Human <- exp(Human.MCMC$pars) %>% apply(2,sd)

dis.Human <- matrix(ncol = 50000, nrow = 17)       # Create an empty matrix
rownames(dis.Human) <- names(M.Human)[1:17]        # The first 17 parameters, one in a row


for (i in 1:50000){
  for (j in 1:17){
    dis.Human[j,i] <- rnorm (1, mean = M.Human[j], sd= SD.Human[j]) # Resample once each parameter for 50000 iterations to get a smooth distribution
  }
}

dis.Human<-melt(dis.Human) # reshape the table, each iteration number now in a colunm, make plotting more convenient
names(dis.Human)=c("Par","Species","Value") # Species is the column name of iteration number
dis.Human$Species = c("Human") # Change the interation number in the Species column to Human, make plotting more convenient
which_Km    <- grep("Km", dis.Human$Par) # The unit of Km in humans is ug/L, which is different from the ug/mL for other species, so a unit conversion is done here.

dis.Human [which_Km,]$Value = dis.Human[which_Km,]$Value/1000

## Rat posteiror distributions
M.Rat  <- exp(Rat.MCMC$pars) %>% apply(2,mean)
SD.Rat <- exp(Rat.MCMC$pars) %>% apply(2,sd)
dis.Rat <- matrix(ncol = 50000, nrow = 17)
rownames(dis.Rat) <- names(M.Rat)[1:17]

for (i in 1:50000){
  for (j in 1:17){
    dis.Rat[j,i] <- rnorm (1, mean = M.Rat[j], sd= SD.Rat[j])
  }
}

dis.Rat<-melt(dis.Rat)
names(dis.Rat)=c("Par","Species","Value")
dis.Rat$Species = c("Rat")

## Mouse posteiror distributions
M.Mouse  <- exp(Mouse.MCMC$pars) %>% apply(2,mean)
SD.Mouse <- exp(Mouse.MCMC$pars) %>% apply(2,sd)
dis.Mouse <- matrix(ncol = 50000, nrow = 17)
rownames(dis.Mouse) <- names(M.Mouse)[1:17]

for (i in 1:50000){
  for (j in 1:17){
    dis.Mouse[j,i] <- rnorm (1, mean = M.Mouse[j], sd= SD.Mouse[j])
  }
 }

dis.Mouse<-melt(dis.Mouse)
names(dis.Mouse)=c("Par","Species","Value")
dis.Mouse$Species = c("Mouse")

## Monkey posteiror distributions
M.Monkey  <- exp(Monkey.MCMC$pars) %>% apply(2,median)
SD.Monkey <- exp(Monkey.MCMC$pars) %>% apply(2,sd)
dis.Monkey <- matrix(ncol = 50000, nrow = 17)
rownames(dis.Monkey) <- names(M.Monkey)[1:17]

for (i in 1:50000){
  for (j in 1:17){
    dis.Monkey[j,i] <- rnorm (1, mean = M.Monkey[j], sd= SD.Monkey[j])
  }
}

dis.Monkey<-melt(dis.Monkey)
names(dis.Monkey)=c("Par","Species","Value")
dis.Monkey$Species = c("Monkey")

## Summary the mouse, rat, monkey, and human posterior distribution in the "Pos.mean" object
Pos.mean <- rbind.data.frame (dis.Mouse,dis.Rat,dis.Monkey,dis.Human)
Pos.mean$log.value <- log(Pos.mean$Value)


## Fig. 3: Densities of posterior parameter uncertainty distributions
## P1: All species gather together in one plot
windowsFonts(Times=windowsFont("Times New Roman")) # Abbreviate the font Times New Roman as Times
Pos.mean$ggjoy = c("A") # Create a new column ggjoy in Pos.mean table. Without this, the plot p1 will have four layers for four species overlapped in each panel, which does not look good.
# With this, y axis is A for all species in one layer in one panel.

Pos.mean$Par = factor(Pos.mean$Par, levels = theta.names[1:17])

p1 = 
  ggplot (Pos.mean, aes(x = as.numeric(log.value), y = as.factor(ggjoy),fill = Species)) + # fill, fill one color for each species
  geom_joy (scale = 8, size = 0.25, rel_min_height = 0.01, alpha = 0.4) + # over size of the entire ggjoy plot
  #geom_joy (data=filter(Pos.mean, Species == "Monkey"),scale = 28, size = 0.25, rel_min_height = 0.01, fill="red", alpha=0.3) +
  #geom_joy (data=filter(Pos.mean, Species == "Mouse"), scale = 28, size = 0.25, rel_min_height = 0.01, fill="green", alpha=0.3) +
  #geom_joy (data=filter(Pos.mean, Species == "Rat"),   scale = 28, size = 0.25, rel_min_height = 0.01, fill="pink", alpha=0.3)+
  scale_fill_manual(values = c("white", "grey", "blue","red")) + 
  scale_colour_manual(values = c("white", "grey", "blue","red"))+ 
  facet_wrap ( ~ Par, nrow = 6, ncol= 3,scale="free_x") # arrange the plot by parameter name, free_x means that the x scale is not fixed, so automatically updated.

p1=
  p1 + theme_joy()+ # change the theme of the plot p1 (incorporate p1 into the theme)
  theme (
    plot.background         = element_blank(),
    text                    = element_text (family = "Times",face="bold",size = 18),
    #panel.border            = element_rect (colour = "black", fill=NA, size=2),
    #panel.background        = element_blank(),
    panel.background        = element_rect (fill   = "#f7f7f7"),
    #panel.grid.major.y      = element_line (size   = 0.5, colour = "grey"),
    #panel.grid.minor.y      = element_blank(),
    panel.grid.major.x      = element_blank(), 
    #panel.grid.minor.x      = element_blank(), 
    #axis.text               = element_text (size   = 18, colour = "black",face="bold"),  
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    axis.ticks.x            = element_line (size   = 1.5, colour = "black"), 
    axis.ticks.y            = element_blank(), 
    axis.text.y             = element_blank(),
    strip.text              = element_text (size   = 18),
    legend.title            = element_text (size   = 18, face="bold"),
    legend.justification    = "none",
    legend.position         = "none",
    legend.text             = element_text (size = 18,face="bold")) + 
    labs (x = "", y = "") 

# Save the image by using the Export button, adjust the height and width, choose image format, preview
ggsave("Fig.3.tiff",scale = 1,
       plot = p1,
       path = "D:/Desktop",
       width = 25, height = 20, units = "cm",dpi=320)


dev.off()






