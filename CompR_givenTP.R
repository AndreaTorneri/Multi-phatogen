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
rho= as.numeric(args[6]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",rho=",rho)
alpha= as.numeric(args[7]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",alpha=",alpha)
pathogen= args[8] # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",pathogen=",pathogen)
ctc.dec= as.numeric(args[9]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",ctc.dec=",ctc.dec)
compl= as.numeric(args[10]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",compl=",compl)
rho.new= as.numeric(args[11]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",rho.new=",rho.new)
alpha.new= as.numeric(args[12]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",alpha.new=",alpha.new)
ctc.dec.new= as.numeric(args[13]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",ctc.dec.new=",ctc.dec.new)
compl.new= as.numeric(args[14]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",compl.new=",compl.new)


if (netw=="ERGM"){
  load("sim_basis_complete_n_1000.RData")
  HH.networks<-HH_sim
  name.s<-paste("TP_ERGM", "_R",R,"_ratioqhqg",ratio.qhqg, "_rho",rho,"_alpha",alpha,"_pathogen",pathogen,"_cdec",ctc.dec,"_comp",compl ,".RData",sep = "")
  load(name.s)
}

if (netw=="Synth"){
  name<-paste("HH_Networks","_nVertex",n.vertex,"_nNetw",n.networks,".RData",sep = "")
  load(name)
  name.s<-paste("TP_Synth_nVertex",n.vertex,"_nNetw",n.networks, "_R",R,"_ratioqhqg",ratio.qhqg , "_rho",rho,"_alpha",alpha,"_pathogen",pathogen,"_cdec",ctc.dec,"_comp",compl,".RData",sep = "")
  load(name.s)
}

lambda.h<-3.34 #average number of daily within household (Mossong et al. 2008 - Belgium)
lambda.g<-8.29 #average number of daily contacts at a community level (Mossong et al. 2008 - Belgium)

transm.prms<-data.frame("beta.g"=inf.path.g*lambda.g, "beta.h"=inf.path.h*lambda.h)

library("network")
source("R_comp_netw.R")


nSim<-100
nSeed<-3082021
set.seed(nSeed)
R.star<-NULL

for (i in 1:nSim){
  temp.HH.netw<-HH.networks[[sample(1:length(HH.networks),1)]]
  R.star<-c(R.star,R0.computation.Inf(HH.network = temp.HH.netw,beta.g = beta.g,nSim = 1,beta.h = beta.h,prob.asym = prob.asym,asymp.rel.inf = asymp.rel.inf,lambda.h = lambda.h,pathogen = pathogen,ctc.dec = ctc.dec,compl = compl))
}

name.r<-paste("RSt_infh",inf.path.h,"_infg",inf.path.g, "_rho",rho.new,"_alpha",alpha.new,"_pathogen",pathogen.new,"_cdec",ctc.dec.new,"_comp",compl.new,".RData",sep = "")

save(R.star, file = name.r)


