################################################################
#### This script computes the transmission parameters for a ####
#### given network. It makes use of the functions defined   ####
#### in the script 'R_comp_netw.R'.                         ####
################################################################

######################################
### Extract Command Line arguments ###
### Set parameter values           ###
######################################

args <- commandArgs(trailingOnly = TRUE) # only arguments after --args are to be returned

# working directory
path_out = args[1] 
cat(",path_out=",path_out) # concatenate and print

# number of cores to run in parallel
cores = as.numeric(args[2]) 
cat(",cores=",cores) # concatenate and print

# type of household network considered - Synthetic or ERGM
netw = as.character(args[3]) 
cat(",netw=",netw)

# number of vertices
n.vertex = as.numeric(args[4])
cat(",n.vertex=",n.vertex)

# number of simulated networks
n.networks = as.numeric(args[5])
cat(",n.networks=",n.networks)

# target value for R
R = as.numeric(args[6])
cat(",R=",R)

# ratio transmission probabilities household contact over global contacts
ratio_tp_hg= as.numeric(args[7])
cat(",ratio_qh_qg=",ratio_qh_qg)

# probability of being symptomatic
rho = as.numeric(args[8]) 
cat(",rho=",rho)

# relative infectiousness of asymptomatic cases
alpha= as.numeric(args[9]) 
cat(",alpha=",alpha)

# the pathogen
pathogen= args[10] 
cat(",pathogen=",pathogen)

# 
ctc.dec= as.numeric(args[11]) 
cat(",ctc.dec=",ctc.dec)

# Behavior change value set to identify transmission rates (household and global) linked to a specific R*
compl= as.numeric(args[12]) 
cat(",compl=",compl)


########################
### TESTING SCENARIO ###
########################

netw = "Synth"
n.vertex = 100
n.networks = 1
R = 3.3
ratio_qh_qg = 8.27
rho = 0.69
alpha = 0.5
pathogen = "COVID-19"
ctc.dec = 1
compl = 0

###################################################
### Load networks, set necessary parameters and ###
### load packages and functions                 ###
###################################################

### Load libraries and functions
library("network")
setwd("C:/Users/LUCP13441/Documents/GitHub/Multi-phatogen")
source("function.TransmParams.R")

### Load network

if(netw == "ERGM"){
  load("sim_basis_complete_n_1000.RData")
  HH.networks <- HH_sim # network data
  name.tp <- paste("TP_ERGM", "_R", R ,"_ratioqhqg", ratio_hg, "_rho", rho,
                   "_alpha", alpha, "_pathogen", pathogen, "_cdec", ctc.dec,
                   "_comp", compl , ".RData", sep = "")
}

if (netw=="Synth"){
  name.network <- paste("HH_Networks", "_nVertex", n.vertex, "_nNetw", n.networks, 
                        ".RData", sep = "")
  setwd("G:/My Drive/PhD/IBM Multipathogen/R-Code")
  load(name.network)
  name.tp <- paste("TP_Synth_nVertex", n.vertex, "_nNetw", n.networks, "_R", R, 
                   "_ratio_qh_qg", ratio_qh_qg, "_rho", rho, "_alpha", alpha, 
                   "_pathogen", pathogen, "_cdec", ctc.dec, "_comp", compl, 
                   ".RData", sep = "")
}



### Set patameters
# Average number of daily contacts
lambda.h <- 3.34 # within household (Mossong et al. 2008 - Belgium)
lambda.g <- 8.29 # at a community level (Mossong et al. 2008 - Belgium)
# the ratio of household contacts to global contacts
ratio_hg <- lambda.h/lambda.g*ratio_qh_qg
# target value
R.rif <- R
# number of simulations
nSim <- 20
# tolerance is 1 % of the target value
tol <- 0.01 * R #R.rif 

###################
### Computation ###
###################

nSeed <- 3082021
set.seed(nSeed)
#trs.prms<-R0.comp.Inf(ratio_hg=ratio_hg, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho),asymp.rel.inf=alpha,lambda.h = lambda.h,pathogen=pathogen,ctc.dec=ctc.dec,compl=compl)
#trs.prms.2<-R0.comp.Inf.bc(ratio_hg=ratio_hg, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho),asymp.rel.inf=alpha,lambda.g = lambda.g,pathogen=pathogen,ctc.dec=ctc.dec,compl=compl)

trs.prms.2 <- R0.comp.RM(ratio_hg = ratio_hg, 
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
inf.path.h<-trs.prms.2$beta.h
inf.path.g<-trs.prms.2$q.g


save(inf.path.h,inf.path.g, file = name.tp)



