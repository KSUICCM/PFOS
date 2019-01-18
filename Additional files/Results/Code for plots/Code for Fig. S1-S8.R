########################################## Fig.S1 ################################# 
# Trace plot using bayesplot for mouse, rat, monkey and human                     #
###################################################################################
## loading R packages
library(magrittr)  # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(reshape)   # melt function to reshape the table
library(readr)     # Used to read RDS data file
library(ggplot2)   # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(bayesplot) # Used to plot the trace plot. Needed for the MCMC objects
library(gridExtra) # for plotting the figure
library(coda)     # for convergence diagnosis

## Loading human, rat, mouse, monkey MCMC combined chains data
monkey.comb  <-readRDS(file="Monkey.comb.rds")
mouse.comb   <-readRDS(file="Mouse.comb.rds")
rat.comb     <-readRDS(file="Rat.comb.rds")
human.comb   <-readRDS(file="Human.comb.rds")
theta.names  <-readRDS(file ="theta.names.rds")

## Estimate the R ratio
R.mouse      <-gelman.diag (mouse.comb)        # Gelman convergence diagnosis
R.rat        <-gelman.diag (rat.comb)          # Gelman convergence diagnosis
R.monkey     <-gelman.diag (monkey.comb)       # Gelman convergence diagnosis
R.human      <-gelman.diag (human.comb)        # Gelman convergence diagnosis

## Convergences plot for Mouse
# Mouse trace plot
p1.mouse <-
  mcmc_trace (
    mouse.comb,
    pars = theta.names[1:17],
    size = 0.5, 
    facet_args = list(ncol= 6, nrow = 3)) + legend_move("none") 

# Mouse probability density function plot
p2.mouse <-
  mcmc_dens_overlay (
    mouse.comb,
    pars = theta.names[1:17],
    size = 0.5,
    facet_args = list(ncol= 6, nrow = 3))+legend_move("none")

grid.arrange(p1.mouse, p2.mouse, nrow = 2)

# Rat trace plot
p1.rat <-
  mcmc_trace (
    rat.comb,
    pars = theta.names[1:17],
    size = 0.5, 
    facet_args = list(ncol= 6, nrow = 3)) + legend_move("none")

# Rat probability density function plot
p2.rat <-
  mcmc_dens_overlay (
    rat.comb,
    pars = theta.names[1:17],
    size = 0.5,
    facet_args = list(ncol= 6, nrow = 3))+legend_move("none")

grid.arrange(p1.rat, p2.rat, nrow = 2)


# Monkey tract plot
p1.monkey <-
  mcmc_trace (
    monkey.comb,
    pars = theta.names[1:17],
    size = 0.5, 
    facet_args = list(ncol= 6, nrow = 3))+ legend_move("none")

# Monkey probability density function plot
p2.monkey <-
    mcmc_dens_overlay (
    monkey.comb,
    pars = theta.names[1:17],
    size = 0.5,
    facet_args = list(ncol= 6, nrow = 3))+legend_move("none")

grid.arrange(p1.monkey, p2.monkey, nrow = 2)

# Human trace plot
p1.human <-
  mcmc_trace (
    human.comb,
    pars = theta.names.h[1:17],
    size = 0.5, 
    facet_args = list(ncol= 6, nrow = 3))+legend_move("none")

# Human probability density function plot
p2.human <-
    mcmc_dens_overlay (
    human.comb,
    pars = theta.names.h[1:17],
    size = 0.5,
    facet_args = list(ncol= 6, nrow = 3))+legend_move("none")

grid.arrange(p1.human, p2.human, nrow = 2)
