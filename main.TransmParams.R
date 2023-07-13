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
R= as.numeric(args[6]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R=",R)
ratio.qhqg= as.numeric(args[7]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",ratio.qhqg=",ratio.qhqg)
rho= as.numeric(args[8]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",rho=",rho)
alpha= as.numeric(args[9]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",alpha=",alpha)
pathogen= args[10] # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",pathogen=",pathogen)
ctc.dec= as.numeric(args[11]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",ctc.dec=",ctc.dec)
compl= as.numeric(args[12]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",compl=",compl)


if (netw=="ERGM"){
  load("sim_basis_complete_n_1000.RData")
  HH.networks<-HH_sim
  name.s<-paste("TP_ERGM", "_R",R,"_ratioqhqg",ratio.qhqg, "_rho",rho,"_alpha",alpha,"_pathogen",pathogen,"_cdec",ctc.dec,"_comp",compl ,".RData",sep = "")
  
}

if (netw=="Synth"){
  name<-paste("HH_Networks","_nVertex",n.vertex,"_nNetw",n.networks,".RData",sep = "")
  load(name)
  name.s<-paste("TP_Synth_nVertex",n.vertex,"_nNetw",n.networks, "_R",R,"_ratioqhqg",ratio.qhqg , "_rho",rho,"_alpha",alpha,"_pathogen",pathogen,"_cdec",ctc.dec,"_comp",compl,".RData",sep = "")
  
}


lambda.h<-3.34 #average number of daily within household (Mossong et al. 2008 - Belgium)
lambda.g<-8.29 #average number of daily contacts at a community level (Mossong et al. 2008 - Belgium)

library("network")

source("R_comp_netw.R")
ratio_hhgl<-lambda.h/lambda.g*ratio.qhqg
R.rif<-R
nSim<-300
tol<-0.005*R.rif #tolerance is 0.5% of the target value
nSeed<-3082021
set.seed(nSeed)
#trs.prms<-R0.comp(ratio_hhgl=ratio.qhqg, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho),asymp.rel.inf=alpha,lambda.h = lambda.h)
trs.prms<-R0.comp.Inf(ratio_hhgl=ratio_hhgl, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho),asymp.rel.inf=alpha,lambda.h = lambda.h,pathogen=pathogen,ctc.dec=ctc.dec,compl=compl)


#load data
inf.path.h<-trs.prms$beta.h/lambda.h
inf.path.g<-trs.prms$beta.g/lambda.g

save(inf.path.h,inf.path.g, file = name.s)



