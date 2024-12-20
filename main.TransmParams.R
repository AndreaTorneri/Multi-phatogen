################################################################
#### This script computes the transmission parameters for a ####
#### given network. It makes use of the functions defined   ####
#### in the script 'R_comp_netw.R'.                         ####
################################################################

### These parameters come from an external file 'TransmParams.csv' 
args <- commandArgs(trailingOnly = TRUE)
# working directory
out <- args[1] 
cat(",out=",out)
# number of cores to run in parallel
cores <- as.numeric(args[2]) 
cat(",cores=",cores)
# type of network, either "Synth" or "ERGM"
netw <- as.character(args[3]) 
cat(",netw=",netw)
# number of individuals in the network
n.vertex <- as.numeric(args[4]) 
cat(",n.vertex=",n.vertex)
# numbers of networks to create
n.networks <- as.numeric(args[5]) #time at which pathogen 2 is insert in the population
cat(",n.networks=",n.networks)
# basic reproductive ratio (?)
R <- as.numeric(args[6]) 
cat(",R=",R)
# ratio between household and global contact rates
ratio.qhqg <- as.numeric(args[7]) 
cat(",ratio.qhqg=",ratio.qhqg)
# long-term interaction (quantity affecting the probability that an individual is infected with 
# a pathogen after experiencing infection with the same or another pathogen)
rho <- as.numeric(args[8]) 
cat(",rho=",rho)
# ?
alpha <- as.numeric(args[9]) 
cat(",alpha=",alpha)
# the pathogen, can be one of the following values:
# "XA", "XS", "XP", "FLU-A", "COVID-19", "RSV", "DELTA", "OMICRON"
pathogen <- args[10] 
cat(",pathogen=",pathogen)
#
ctc.dec <- as.numeric(args[11]) 
cat(",ctc.dec=",ctc.dec)
# 
compl <- as.numeric(args[12]) 
cat(",compl=",compl)

### For testing purposes, use this set of parameters:
netw <- "Synth"
n.vertex <- 100
n.networks <- 1
R <- 1.5
ratio.qhqg <- 1
rho <- 1
alpha <- 1
pathogen <- "XS"
ctc.dec <- 1
compl <- 1

### Load the corresponding networks
if (netw == "ERGM"){
  load("sim_basis_complete_n_1000.RData")
  HH.networks <- HH_sim
  name.s <- paste("TP_ERGM", "_R",R,"_ratioqhqg", ratio.qhqg, "_rho",rho, "_alpha", alpha, 
                  "_pathogen", pathogen, "_cdec", ctc.dec, "_comp", compl , ".RData", sep = "")
}
if (netw == "Synth"){
  name <- paste("HH_Networks", "_nVertex", n.vertex, "_nNetw", n.networks, ".RData", sep = "")
  setwd("C:/Users/LUCP13441/Documents/GitHub/Multi-phatogen")
  load(name)
  name.s <- paste("TP_Synth_nVertex", n.vertex, "_nNetw", n.networks, "_R", R, "_ratioqhqg", ratio.qhqg,
                  "_rho", rho, "_alpha", alpha, "_pathogen", pathogen, "_cdec", ctc.dec, "_comp",
                  compl, ".RData", sep = "")
}

# Average number of daily within household (Mossong et al. 2008 - Belgium)
lambda.h <- 3.34 
# Average number of daily contacts at a community level (Mossong et al. 2008 - Belgium)
lambda.g <- 8.29 

# load the necessary files and packages
library("network")
source("R_comp_netw.R")


ratio_hhgl <- lambda.h/lambda.g*ratio.qhqg
R.rif <- R
nSim <- 20
tol <- 0.01*R.rif # tolerance is 0.5% of the target value
nSeed <- 3082021
set.seed(nSeed)
#trs.prms<-R0.comp.Inf(ratio_hhgl=ratio_hhgl, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho),asymp.rel.inf=alpha,lambda.h = lambda.h,pathogen=pathogen,ctc.dec=ctc.dec,compl=compl)
#trs.prms.2<-R0.comp.Inf.bc(ratio_hhgl=ratio.qhqg, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho),asymp.rel.inf=alpha,lambda.g = lambda.g,pathogen=pathogen,ctc.dec=ctc.dec,compl=compl)
trs.prms.2 <- R0.comp.RM(ratio_hhgl = ratio.qhqg, 
                         HH.network = HH.networks, 
                         nSim = nSim, 
                         tol = tol,
                         R.rif = R.rif, 
                         prob.asym = (1-rho),
                         asymp.rel.inf = alpha,
                         lambda.g = lambda.g,
                         pathogen = pathogen,
                         ctc.dec = ctc.dec,
                         compl = compl)


#load data
inf.path.h <- trs.prms.2$beta.h
inf.path.g <- trs.prms.2$q.g


save(inf.path.h,inf.path.g, file = name.s)



