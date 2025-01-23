# Input parameters from csv
args <- commandArgs(trailingOnly = TRUE)
# working directory
out = args[1] 
cat(",out=",out)
# number of cores to run in parallel
cores = as.numeric(args[2]) 
cat(",cores=",cores)
# time at which pathogen 2 is introduced in the population
t2 = as.numeric(args[3]) 
cat(",t2=",t2)
# short-term interaction parameter: acquiring 2 while having 1 (if >1 cooperative effect - if <1 competing)
sigma12 = as.numeric(args[4]) 
cat(",sigma12=",sigma12)
# short-term interaction parameter: acquiring 1 while having 2 (if >1 cooperative effect - if <1 competing)
sigma21 = as.numeric(args[5]) 
cat(",sigma21=",sigma21)
# proportion of immune cases (not used at the moment)
prop.immune = as.numeric(args[6]) 
cat(",prop.immune=",prop.immune)
# number of initial cases for path 1
nSeeds.1= as.numeric(args[7]) 
cat(",nSeeds.1=",nSeeds.1)
# number of initial cases for path 2
nSeeds.2= as.numeric(args[8]) 
cat(",nSeeds.2=",nSeeds.2)
# probability of being symptomatic for path 1
rho.1= as.numeric(args[9]) 
cat(",rho.1=",rho.1)
# probability of being symptomatic for path 2
rho.2= as.numeric(args[10]) 
cat(",rho.2=",rho.2)
# relative infectiousness of asymptomatic cases (pathogen1)
alpha.as.1= as.numeric(args[11]) 
cat(",alpha.as.1=",alpha.as.1)
# relative infectiousness of asymptomatic cases (pathogen2)
alpha.as.2= as.numeric(args[12]) 
cat(",alpha.as.2=",alpha.as.2)
# type of household network considered - Synthetic or ERGM
netw = as.character(args[13]) 
cat(",netw=",netw)
# number of vertexes 
n.vertex = as.numeric(args[14]) 
cat(",n.vertex=",n.vertex)
# number of simulated networks
n.networks = as.numeric(args[15]) 
cat(",n.networks=",n.networks)
# Reproduction number path 1 (household R*)
R.1= as.numeric(args[16]) 
cat(",R.1=",R.1)
# Reproduction number path 2 (household R*)
R.2= as.numeric(args[17]) 
cat(",R.2=",R.2)
# ratio transmission probability given household contact over  global contacts
ratio.qhqg= as.numeric(args[18]) 
cat(",ratio.qhqg=",ratio.qhqg)
# long-term interaction parameter: acquiring 2 while having experienced (and recovered from) 1 
lli.1= as.numeric(args[19]) 
cat(",lli.1=",lli.1)
# long-term interaction parameter: acquiring 1 while having experienced (and recovered from) 2
lli.2= as.numeric(args[20])  
cat(",lli.2=",lli.2)
# character variable identifying pathogen 1
pathogen.1= as.character(args[21]) 
cat(",pathogen.1=",pathogen.1)
# character variable identifying pathogen 2
pathogen.2= as.character(args[22]) 
cat(",pathogen.2=",pathogen.2)
# parameter multiplying the household contact rate after home isolation
contact.reduction= as.numeric(args[23]) 
cat(",contact.reduction=",contact.reduction)
# time at which simulations stop
t.stop= as.numeric(args[24]) 
cat(",t.stop=",t.stop)
# time of additional seeding
t.seed= as.numeric(args[25]) 
cat(",t.seed=",t.seed)
# proportion of individuals changing behavior (home isolation) after being infected with pathogen 1
bc.1= as.numeric(args[26]) 
cat(",bc.1=",bc.1)
# proportion of individuals changing behavior (home isolation) after being infected with pathogen 2
bc.2= as.numeric(args[27]) 
cat(",bc.2=",bc.2)
# Boolean identifying whether someone can be re-infected with the same pathogen (1 yes, 0 no)
reinf= as.numeric(args[28]) 
cat(",reinf=",reinf)
# ID for different type of waning of immunity
typeIC= as.numeric(args[29]) 
cat(",typeIC=",typeIC)
# contact reduction value set to identify transmission rates (household and global) linked to a specific R*
contact.reduction.TP= as.numeric(args[30]) 
cat(",contact.reduction.TP=",contact.reduction.TP)
# behavior change value (for pathogen 1) set to identify transmission rates (household and global) linked to a specific R*
bc.1.TP= as.numeric(args[31]) 
cat(",bc.1.TP=",bc.1.TP)
# behavior change value (for pathogen 2) set to identify transmission rates (household and global) linked to a specific R*
bc.2.TP= as.numeric(args[32]) 
cat(",bc.2.TP=",bc.2.TP)
# Boolean for heterologous effects (1 yes 0 no) - Not used currently
het.vac= as.numeric(args[33]) 
cat(",het.vac=",het.vac)
# parameter to define the length of immunity that have the same overall "effect" (area underneath the curve)
t.imm.lim= as.numeric(args[34]) 
cat(",t.imm.lim=",t.imm.lim)
# Decrease in the  number of global contact rates compared to baseline
dec.gc=as.numeric(args[35]) 
cat(",dec.gc=",dec.gc)

### For testing purposes, use this set of parameters:
t2 = 0 # time at which pathogen 2 is introduced in the population          
sigma12 = 1 # short-term interaction parameter: acquiring 2 while having 1 (if >1 cooperative effect - if <1 competing)                      
sigma21 = 1 # short-term interaction parameter: acquiring 1 while having 2 (if >1 cooperative effect - if <1 competing) 
prop.immune = 0 # proportion of immune cases (not used at the moment)
nSeeds.1 = 20 # number of initial cases for path 1
nSeeds.2 = 20 # number of initial cases for path 2
rho.1 = 0.69 # probability of being symptomatic for path 1
rho.2 = 0.67 # probability of being symptomatic for path 2
alpha.as.1 = 0.5 # relative infectiousness of asymptomatic cases (pathogen1)
alpha.as.2 = 0.33 # relative infectiousness of asymptomatic cases (pathogen2)
netw = "Synth" # type of household network considered - Synthetic or ERGM
n.vertex = 100 # number of vertices 
n.networks = 1 # number of simulated networks
R.1 = 3.3 # reproduction number path 1 (household R*)
R.2 = 1.3 # reproduction number path 2 (household R*)
ratio.qhqg = 8.27 # ratio transmission probability given household contact over global contacts
lli.1 = 1 # long-term interaction parameter: acquiring 2 while having experienced (and recovered from) 1 
lli.2 = 1 # long-term interaction parameter: acquiring 1 while having experienced (and recovered from) 2
pathogen.1 = "COVID-19" # character variable identifying pathogen 1
pathogen.2 = "FLU-A" # character variable identifying pathogen 2
contact.reduction = 1 # parameter multiplying the household contact rate after home isolation
t.stop = 365 # time at which simulations stop
t.seed = 1000 # time of additional seeding
bc.1 = 0 # proportion of individuals changing behavior (home isolation) after being infected with pathogen 1
bc.2 = 0 # proportion of individuals changing behavior (home isolation) after being infected with pathogen 2
reinf = 0 # boolean identifying whether someone can be re-infected with the same pathogen (1 yes, 0 no)
typeIC = 1  # ID for different type of waning of immunity
contact.reduction.TP = 1 # contact reduction value set to identify transmission rates (household and global) linked to a specific R*
bc.1.TP = 0 # behavior change value (for pathogen 1) set to identify transmission rates (household and global) linked to a specific R*
bc.2.TP = 0 # behavior change value (for pathogen 2) set to identify transmission rates (household and global) linked to a specific R*
het.vac = 1 # boolean for heterologous effects (1 yes 0 no) - Not used currently
t.imm.lim = 10 # parameter to define the length of immunity that have the same overall "effect" (area underneath the curve)
dec.gc = 1 # decrease in the  number of global contact rates compared to baseline

### Load necessary packages
library(ergm)
library(RGeode)

#Two type of household networks can be loaded (ERGM (data-driven) - Synthetic (household size representative but random mixing))
#For now considered only Synthetic networks
# To note, the type of network together with other characteristics (e.g., R*), will give you the value of the transmission parameters for global and local contacts
# this values can be computed with another Rscript present in the repo - mainTransmParams.

if (netw=="ERGM"){
  load("sim_basis_complete_n_1000.RData")
  HH.networks<-HH_sim
  name.s<-paste("TransParam_ERGMNetworks_nVertex",n.vertex,"_nNetw",n.networks, "_R",R.1,"_ratioqhqg",ratio.qhqg, "_rho",rho.1,"_alpha",alpha.as.1,".RData",sep = "")
  load(name.s)
  inf.path.1.h<-inf.path.h
  inf.path.1.g<-inf.path.g
  name.s<-paste("TransParam_ERGMNetworks_nVertex","_nNetw",n.networks, "_R",R.2,"_ratioqhqg",ratio.qhqg, "_rho",rho.2,"_alpha",alpha.as.2,".RData",sep = "")
  load(name.s)
  inf.path.2.h<-inf.path.h
  inf.path.2.g<-inf.path.g
}

if (netw=="Synth"){
  # load the corresponding network
  name.n <- paste("HH_Networks","_nVertex",n.vertex,"_nNetw",n.networks,".RData",sep = "")
  setwd("G:/My Drive/PhD/IBM Multipathogen/R-Code")
  load(name.n)
  
  # set pathogen.1
  pathogen.TP <- pathogen.1
  if(pathogen.1 == "DELTA" | pathogen.1=="OMICRON"){
    pathogen.TP <- "COVID-19"
  }
  if(pathogen.1 == "XP" | pathogen.1 == "XA"){
    pathogen.TP <- "XS" 
  }
  
  # load the corresponding 
  name.s <- paste("TP_Synth_nVertex",n.vertex,"_nNetw",n.networks, "_R",R.1,"_ratioqhqg",ratio.qhqg , "_rho",rho.1,"_alpha",alpha.as.1,"_pathogen",pathogen.TP,"_cdec",contact.reduction.TP,"_comp",bc.1.TP,".RData",sep = "")
  load(name.s)
  inf.path.1.h<-inf.path.h
  inf.path.1.g<-inf.path.g
  pathogen.TP<-pathogen.2
  if (pathogen.2=="DELTA" | pathogen.2=="OMICRON"){ pathogen.TP<-"COVID-19" }
  if (pathogen.2=="XP" | pathogen.2=="XA"){ pathogen.TP<-"XS" }
  name.s<-paste("TP_Synth_nVertex",n.vertex,"_nNetw",n.networks, "_R",R.2,"_ratioqhqg",ratio.qhqg , "_rho",rho.2,"_alpha",alpha.as.2,"_pathogen",pathogen.TP,"_cdec",contact.reduction.TP,"_comp",bc.2.TP,".RData",sep = "")
  load(name.s)
  inf.path.2.h<-inf.path.h
  inf.path.2.g<-inf.path.g
}

#Mean number of daily contact at a global level (Using SOCRATES 15/02/2022 )
lambda.g<-8.29*dec.gc


#Compute the reproduction number related to the selected network. 
source("function.multipathogen.new.R")
nSim<-100
epi.outbreak<-list()
nSeed<-1062021
set.seed(nSeed)

nm<-paste("t2_",t2, "_sigma12_",sigma12,"_sigma21_",sigma21,"_qh1_",inf.path.1.h,"_qg1_",inf.path.1.g,"_qh2_",inf.path.2.h,"_qg2_",inf.path.2.g, "_rho1_",rho.1,"_rho2_",rho.2,"_alpha1_",alpha.as.1,"_alpha2_",alpha.as.2,"_Path1",pathogen.1,"_Path2",pathogen.2,"_lli.1",lli.1,"_lli2",lli.2, sep = "")
print(nm)
for (i in 1:nSim){
  print(i)
  temp.HH.netw<-HH.networks[[sample(1:length(HH.networks),1)]]
  epi.outbreak[[i]]<-sim.multipathogen(HH.network = temp.HH.netw, t2=t2, lambda.g = lambda.g, sigma12 = sigma12, sigma21 = sigma21, prop.immune = prop.immune, nSeeds.1 = nSeeds.1, nSeeds.2 = nSeeds.2, rho.1 = rho.1, rho.2 = rho.2, inf.path.1.h = inf.path.1.h,inf.path.1.g = inf.path.1.g, inf.path.2.h = inf.path.2.h,inf.path.2.g = inf.path.2.g, alpha.as.1=alpha.as.1,alpha.as.2=alpha.as.2, lli.1=lli.1,lli.2=lli.2, pathogen.1=pathogen.1, pathogen.2=pathogen.2, contact.reduction=contact.reduction, t.stop=t.stop, t.seed=t.seed, bc.1=bc.1, bc.2=bc.2, reinf=reinf, typeIC=typeIC, het.vac=het.vac, t.imm.lim=t.imm.lim)
}

scen<-paste(netw,"_nVertex",n.vertex,"_nNetw",n.networks,pathogen.1,"_&_",pathogen.2,sep ="")

name<-paste("MP_",scen,"_R1",R.1,"_R2",R.2,"_qhqg",ratio.qhqg, "_t2",t2,"_sigma12_",sigma12,"_sigma21_",sigma21,"_alpha1",alpha.as.1,"_alpha2",alpha.as.2,"_rho1",rho.1,"_rho2",rho.2,"_lli1",lli.1,"_lli2",lli.2,"_Net",netw,"_CtcRed",contact.reduction,"_PImm",prop.immune,"_tSeed",t.seed, "_bc1",bc.1,"_bc2",bc.2,".RData", sep = "")
setwd(out)
save(epi.outbreak, file = name)

