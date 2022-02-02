args <- commandArgs(trailingOnly = TRUE)
out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
netw = as.character(args[3]) #time at which pathogen 2 is insert in the population
cat(",netw=",netw)
n.vertex = as.numeric(args[4]) #time at which pathogen 2 is insert in the population
cat(",n.vertex=",n.vertex)
n.networks = as.numeric(args[5]) #time at which pathogen 2 is insert in the population
cat(",n.networks=",n.networks)
R.1= as.numeric(args[6]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R.1=",R.1)
R.2= as.numeric(args[7]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R.2=",R.2)
ratio.qhqg= as.numeric(args[8]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",ratio_qhqg=",ratio_qhqg)

ratio_hhgl

if (netw=="ERGM"){
  load("sim_basis_complete_n_1000.RData")
  HH.networks<-HH_sim
  name.s<-paste("TransParam_ERGMNetworks", "_R1",R.1,"_R2",R.2,"_ratioqhqg",ratio.qhqg ,".RData",sep = "")
  
}
if (netw=="Synth"){
  name<-paste("HH_Networks","_nVertex",n.vertex,"_nNetw",n.networks,".RData",sep = "")
  load(name)
  name.s<-paste("TransParam_SynthNetworks","_nVertex",n.vertex,"_nNetw",n.networks, "_R1",R.1,"_R2",R.2,"_ratioqhqg",ratio.qhqg ,".RData",sep = "")
  
}

lambda.h<-3.34 #average number of daily within household (Mossong et al. 2008 - Belgium)
lambda.g<-8.29 #average number of daily contacts at a community level (Mossong et al. 2008 - Belgium)

source("R_comp_netw.R")
ratio_hhgl<-lambda.h/lambda.g*ratio.qhqg
R.rif<-R.1
nSim<-100
tol<-0.05
nSeed<-3082021
set.seed(nSeed)
trs.prms.path<-R0.comp(ratio_hhgl=ratio_hhgl, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif )

#load data
inf.path.1.h<-trs.prms$beta.h/lambda.h
inf.path.1.g<-trs.prms$beta.g/lambda.g

R.rif<-R.2
trs.prms<-R0.comp(ratio_hhgl=ratio_hhgl, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif )

inf.path.2.h<-trs.prms$beta.h/lambda.h
inf.path.2.g<-trs.prms$beta.g/lambda.g


save(inf.path.1.h,inf.path.1.g,inf.path.2.h,inf.path.2.g, file = name.s)



