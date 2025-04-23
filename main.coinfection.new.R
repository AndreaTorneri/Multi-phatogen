#########################################################
#### This script runs the main simulation of the IBM ####
#########################################################


######################################
### Extract Command Line arguments ###
### Set parameter values           ###
######################################

# Input parameters from csv
args <- commandArgs(trailingOnly = TRUE)

# working directory
out <<- args[1] 
cat(",out=", out)

# number of cores to run in parallel
cores <<- as.numeric(args[2]) 
cat(",cores=", cores)

# time at which pathogen 2 is introduced in the population
t2 <<- as.numeric(args[3]) 
cat(",t2=", t2)

# short-term interaction parameter: acquiring the other disease while having the one (if >1 cooperative effect - if <1 competing)
sigma <<- NULL
sigma[1] <<- as.numeric(args[4]) 
cat(",sigma12=", sigma[1])
sigma[2] <<- as.numeric(args[5]) 
cat(",sigma21=", sigma[2])

# proportion of immune cases (not used at the moment)
prop.immune <<- as.numeric(args[6]) 
cat(",prop.immune=", prop.immune)

# number of initial cases
n.seeds <<- NULL
n.seeds[1] <<- as.numeric(args[7]) 
cat(",nSeeds.1=", n.seeds[1])
n.seeds[2] <<- as.numeric(args[8]) 
cat(",nSeeds.2=", n.seeds[2])

# probability of being symptomatic
rho = NULL 
rho[1] = as.numeric(args[9]) 
cat(",rho.1=", rho[1])
rho[2] = as.numeric(args[10]) 
cat(",rho.2=", rho[2])

# relative infectiousness of asymptomatic cases
alpha.as = NULL
alpha.as[1]= as.numeric(args[11]) 
cat(",alpha.as.1=", alpha.as[1])
alpha.as[2] = as.numeric(args[12]) 
cat(",alpha.as.2=", alpha.as[2])

# type of household network considered - Synthetic or ERGM
netw = as.character(args[13]) 
cat(",netw=", netw)

# number of vertices 
n.vertex = as.numeric(args[14]) 
cat(",n.vertex=", n.vertex)

# number of simulated networks
n.networks = as.numeric(args[15]) 
cat(",n.networks=", n.networks)

# Reproduction number (household R*)
R = NULL
R[1] = as.numeric(args[16]) 
cat(",R.1=", R[1])
R[2] = as.numeric(args[17]) 
cat(",R.2=", R[2])

# ratio transmission probability given household contact over  global contacts
ratio.qhqg = as.numeric(args[18]) 
cat(",ratio.qhqg=", ratio.qhqg)

# long-term interaction parameter: acquiring the other pathogen while having experienced one pathogen
long.int[1] = as.numeric(args[19]) 
cat(",lli.1=", long.int[1])
long.int[2] = as.numeric(args[20])  
cat(",lli.2=", long.int[2])

# character variable identifying
pathogen = NULL
pathogen[1] = as.character(args[21]) 
cat(",pathogen1=", pathogen[1])
pathogen[2] = as.character(args[22]) 
cat(",pathogen2=", pathogen[2])

# parameter multiplying the household contact rate after home isolation
contact.reduction = as.numeric(args[23]) 
cat(",contact.reduction=", contact.reduction)

# time at which simulations stop
t.stop = as.numeric(args[24]) 
cat(",t.stop=", t.stop)

# time of additional seeding
t.seed = as.numeric(args[25]) 
cat(",t.seed=", t.seed)

# proportion of individuals changing behavior (home isolation) after being infected 
behavior.change = NULL
behavior.change[1] = as.numeric(args[26]) 
cat(",bc.1=", behavior.change[1])
behavior.change.[2] = as.numeric(args[27]) 
cat(",bc.2=", behavior.change[2])

# Boolean identifying whether someone can be re-infected with the same pathogen (1 yes, 0 no)
reinfection = as.numeric(args[28]) 
cat(",reinf=", reinfection)

# ID for different type of waning of immunity
typeIC = as.numeric(args[29]) 
cat(",typeIC=", typeIC)

# contact reduction value set to identify transmission rates (household and global) linked to a specific R*
contact.reduction.TP = as.numeric(args[30]) 
cat(",contact.reduction.TP=", contact.reduction.TP)

# behavior change value (for pathogen 1) set to identify transmission rates (household and global) linked to a specific R*
behavior.change.TP = NULL
behavior.change.TP[1] = as.numeric(args[31]) 
cat(",bc.1.TP=", behavior.change.TP[1])
behavior.change.TP[2] = as.numeric(args[32]) 
cat(",bc.2.TP=", behavior.change.TP[2])

# Boolean for heterologous effects (1 yes 0 no) - Not used currently
het.vac = as.numeric(args[33]) 
cat(",het.vac=", het.vac)

# parameter to define the length of immunity that have the same overall "effect" (area underneath the curve)
t.imm.lim = as.numeric(args[34]) 
cat(",t.imm.lim=", t.imm.lim)

# Decrease in global contact rates compared to baseline
decrease.gc =as.numeric(args[35]) 
cat(",dec.gc=", decrease.gc)

vaccine.effectiveness <- as.numeric(args[36]) # TO DO
cat(",vaccine.effectiveness=", vaccine.effectiveness)



########################
### TESTING SCENARIO ###
########################

t2 <- 0 # time at which pathogen 2 is introduced in the population          
sigma <- c(1, 1) # short-term interaction parameter: acquiring 2 while having 1 (if >1 cooperative effect - if <1 competing) 
prop.immune <- 0 # proportion of immune cases (not used at the moment)
n.seeds <- c(2, 2) # number of initial cases for path 1 and 2
rho <- c(0.69, 0.67) # probability of being symptomatic
alpha.as <- c(0.5, 0.33) # relative infectiousness of asymptomatic cases
netw <- "Synth" # type of household network considered - Synthetic or ERGM
n.vertex <- 100 # number of vertices 
n.networks <- 1 # number of simulated networks
R <- c(3.3, 1.3) # reproduction number
ratio.qhqg <- 8.27 # ratio transmission probability given household contact over global contacts
long.int <- c(1, 1) # long-term interaction parameter: acquiring 2 while having experienced (and recovered from) 1
pathogen <- c("COVID-19", "FLU-A")
contact.reduction <- 1 # parameter multiplying the household contact rate after home isolation
t.stop <- 365 # time at which simulations stop
t.seed <- 1000 # time of additional seeding
behavior.change <- c(0.5, 0.25) # proportion of individuals changing behavior (home isolation) after being infected
reinfection <- 0 # boolean identifying whether someone can be re-infected with the same pathogen (1 yes, 0 no)
typeIC <- 0  # ID for different type of waning of immunity
contact.reduction.TP <- 1 # contact reduction value set to identify transmission rates (household and global) linked to a specific R*
behavior.change.TP <- c(0, 0) # behavior change value set to identify transmission rates (household and global) linked to a specific R*
het.vac <- 1 # boolean for heterologous effects (1 yes 0 no) - Not used currently
t.imm.lim <- 10 # parameter to define the length of immunity that have the same overall "effect" (area underneath the curve)
decrease.gc <- 1 # decrease in the  number of global contact rates compared to baseline
vaccination.coverage <- c(0.8, 0.8) # the probability of being vaccinated
prop.vaccinated <- c(0, 0) # the proportion of the population vaccinated at the start of the simulation

###################################################
### Load networks, set necessary parameters and ###
### load packages and functions                 ###
###################################################

library(ergm)
library(RGeode)
library(dplyr)
library(tidyverse)
library(data.table)

# Two type of household networks can be loaded (ERGM (data-driven) - Synthetic (household size representative but random mixing))
# For now considered only Synthetic networks
# To note, the type of network together with other characteristics (e.g., R*), will give you the value of the transmission parameters for global and local contacts
# this values can be computed with another Rscript present in the repo - mainTransmParams.

if(netw == "ERGM"){
  # load the corresponding network
  load("sim_basis_complete_n_1000.RData")
  HH.networks <- HH_sim
  # load and store transmission parameters for each pathogen 
  inf.h <- inf.g <- NULL
  for(p in 1:length(pathogen)){
    name.s <- paste("TransParam_ERGMNetworks_nVertex", n.vertex, "_nNetw", n.networks, "_R", R[p],
                    "_ratioqhqg", ratio.qhqg, "_rho", rho[1], "_alpha", alpha.as[p], ".RData", sep = "")
    load(name.s)
    inf.h[p] <- inf.path.h
    inf.g[p] <- inf.path.g
  }
}

if(netw == "Synth"){
  # load the corresponding network
  name.n <- paste("HH_Networks", "_nVertex", n.vertex, "_nNetw", n.networks, ".RData", sep = "")
  setwd("G:/My Drive/PhD/IBM Multipathogen/R-Code") #only necessary for testing purposes if working on local machine. 
  load(name.n)
  
  # load and store transmission parameters for each pathogen 
  inf.h <- inf.g <- NULL
  for(p in 1:length(pathogen)){
    pathogen.TP <- pathogen[p]
    if(pathogen[p] == "DELTA" | pathogen[p] == "OMICRON"){pathogen.TP <- "COVID-19"}
    if(pathogen[p] == "XP" | pathogen[p] == "XA"){pathogen.TP <- "XS"}
    name.s <- paste("TP_Synth_nVertex", n.vertex, "_nNetw", n.networks, "_R", R[p], "_ratioqhqg", ratio.qhqg,
                    "_rho", rho[p], "_alpha", alpha.as[p], "_pathogen", pathogen.TP, "_cdec", contact.reduction.TP,
                    "_comp", behavior.change.TP[p], ".RData", sep = "")
    load(name.s)
    inf.h[p] <- inf.path.h
    inf.g[p] <- inf.path.g
  }
}

# Mean number of daily contacts at a global level (Using SOCRATES 15/02/2022)
lambda.g <- 8.29 * decrease.gc


# Compute the reproduction number related to the selected network. 
source("C:/Users/LUCP13441/Documents/GitHub/Multi-phatogen/function.multipathogen.new.R")
n.sim <- 3
epi.outbreak <- list()
n.seed <- 1062021
set.seed(n.seed)

#nm <- paste0("t2_", t2, "_sigma12_", sigma[1], "_sigma21_", sigma[2], "_qh1_", inf.h[1], "_qg1_", inf.g[1], "_qh2_",inf.h[2], "_qg2_", inf.g[2], "_rho1_", rho[1], "_rho2_", rho[2], "_alpha1_", alpha.as[1], "_alpha2_", alpha.as[2], "_Path1_", pathogen[1], "_Path2_", pathogen[2], "_lli.1_", long.int[1], "_lli2_", long.int[2])
nm <- paste0("t2: ", t2, ", sigma12: ", sigma[1], ", sigma21: ", sigma[2], ", qh1: ", inf.h[1], ", qg1: ", inf.g[1],
             ", qh2: ",inf.h[2], ", qg2: ", inf.g[2], ", rho1: ", rho[1], ", rho2: ", rho[2], ", alpha1: ", alpha.as[1], 
             ", alpha2: ", alpha.as[2], ", Path1: ", pathogen[1], ", Path2: ", pathogen[2], ", lli.1: ", long.int[1], ", lli2: ", long.int[2])
print(nm)


######################
### SIMULATION IBM ###
######################

for(i in 1:n.sim){
  print(paste0("simulation ", i))
  # select a network at random
  temp.HH.netw <- HH.networks[[sample(1:length(HH.networks), 1)]]
  epi.outbreak[[i]] <- sim.multipathogen(HH.network = temp.HH.netw, t2 = t2, t.seed = t.seed,
                                         lambda.g = lambda.g)#, t2 = t2, lambda.g = lambda.g, 
                                         #prop.immune = prop.immune, sigma = sigma,
                                         #n.seeds = n.seeds, rho = rho, inf.path.h = inf.h, inf.path.g = inf.g,
                                         #alpha.as = alpha.as, long.int = long.int, pathogen = pathogen, 
                                         #contact.reduction = contact.reduction, t.stop = t.stop, 
                                         #t.seed = t.seed, behavior.change = behavior.change,reinfection=reinfection, 
                                         #typeIC = typeIC, het.vac = het.vac, t.imm.lim = t.imm.lim)
}

scen <- paste(netw,"_nVertex",n.vertex,"_nNetw",n.networks,pathogen[1],"_&_",pathogen[2],sep ="")

name <- paste("MP_",scen,"_R1",R.1,"_R2",R.2,"_qhqg",ratio.qhqg, "_t2",t2,"_sigma12_",sigma[1],"_sigma21_",sigma[2],"_alpha1",alpha.as[1],"_alpha2",alpha.as[2],"_rho1",rho[1],"_rho2",rho[2],"_lli1",long.int[1],"_lli2",long.int[2],"_Net",netw,"_CtcRed",contact.reduction,"_PImm",prop.immune,"_tSeed",t.seed, "_bc1",behavior.change.1,"_bc2",behavior.change.2,".RData", sep = "")
setwd(out)
save(epi.outbreak, file = name)

