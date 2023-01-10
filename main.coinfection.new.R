#Input parameters from csv
args <- commandArgs(trailingOnly = TRUE)
out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
t2 = as.numeric(args[3]) #time at which pathogen 2 is insert in the population
cat(",t2=",t2)
sigma12 = as.numeric(args[4]) # short-term interaction parameter: acquiring 2 while having 1 (this value can)
cat(",sigma12=",sigma12)
sigma21 = as.numeric(args[5]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",sigma21=",sigma21)
prop.immune = as.numeric(args[6]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",prop.immune=",prop.immune)
nSeeds.1= as.numeric(args[7]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",nSeeds.1=",nSeeds.1)
nSeeds.2= as.numeric(args[8]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",nSeeds.2=",nSeeds.2)
rho.1= as.numeric(args[9]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",rho.1=",rho.1)
rho.2= as.numeric(args[10]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",rho.2=",rho.2)
alpha.as.1= as.numeric(args[11]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",alpha.as.1=",alpha.as.1)
alpha.as.2= as.numeric(args[12]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",alpha.as.2=",alpha.as.2)
netw = as.character(args[13]) #time at which pathogen 2 is insert in the population
cat(",netw=",netw)
n.vertex = as.numeric(args[14]) #time at which pathogen 2 is insert in the population
cat(",n.vertex=",n.vertex)
n.networks = as.numeric(args[15]) #time at which pathogen 2 is insert in the population
cat(",n.networks=",n.networks)
R.1= as.numeric(args[16]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R.1=",R.1)
R.2= as.numeric(args[17]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R.2=",R.2)
ratio.qhqg= as.numeric(args[18]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",ratio.qhqg=",ratio.qhqg)
lli.1= as.numeric(args[19]) # long-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",lli.1=",lli.1)
lli.2= as.numeric(args[20]) # long-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",lli.2=",lli.2)
pathogen.1= as.character(args[21]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",pathogen.1=",pathogen.1)
pathogen.2= as.character(args[22]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",pathogen.2=",pathogen.2)
contact.reduction= as.numeric(args[23]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",contact.reduction=",contact.reduction)
t.stop= as.numeric(args[24]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",t.stop=",t.stop)
t.seed= as.numeric(args[25]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",t.seed=",t.seed)
bc.1= as.numeric(args[26]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",bc.1=",bc.1)
bc.2= as.numeric(args[27]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",bc.2=",bc.2)
reinf= as.numeric(args[28]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",reinf=",reinf)
typeIC= as.numeric(args[29]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",typeIC=",typeIC)
contact.reduction.TP= as.numeric(args[30]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",contact.reduction.TP=",contact.reduction.TP)
bc.1.TP= as.numeric(args[31]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",bc.1.TP=",bc.1.TP)
bc.2.TP= as.numeric(args[32]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",bc.2.TP=",bc.2.TP)



#Input parameters - fixed
library(ergm)
library(RGeode)


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
  name.n<-paste("HH_Networks","_nVertex",n.vertex,"_nNetw",n.networks,".RData",sep = "")
  load(name.n)
  name.s<-paste("TP_Synth_nVertex",n.vertex,"_nNetw",n.networks, "_R",R.1,"_ratioqhqg",ratio.qhqg , "_rho",rho,"_alpha",alpha,"_pathogen",pathogen,"_cdec",contact.reduction.TP,"_comp",bc.1.TP,".RData",sep = "")
  load(name.s)
  inf.path.1.h<-inf.path.h
  inf.path.1.g<-inf.path.g
  name.s<-paste("TP_Synth_nVertex",n.vertex,"_nNetw",n.networks, "_R",R.2,"_ratioqhqg",ratio.qhqg , "_rho",rho,"_alpha",alpha,"_pathogen",pathogen,"_cdec",contact.reduction.TP,"_comp",bc.2.TP,".RData",sep = "")
  load(name.s)
  inf.path.2.h<-inf.path.h
  inf.path.2.g<-inf.path.g
}

#Mean number of daily contact at a global level (Using SOCRATES 15/02/2022 )
lambda.g<-8.29


#Compute the reproduction number related to the selected network. 
source("function.multipathogen.new.R")
nSim<-50
epi.outbreak<-list()
nSeed<-1062021
set.seed(nSeed)

nm<-paste("t2_",t2, "_sigma12_",sigma12,"_sigma21_",sigma21,"_qh1_",inf.path.1.h,"_qg1_",inf.path.1.g,"_qh2_",inf.path.2.h,"_qg2_",inf.path.2.g, "_rho1_",rho.1,"_rho2_",rho.2,"_alpha1_",alpha.as.1,"_alpha2_",alpha.as.2,"_Path1",pathogen.1,"_Path2",pathogen.2,"_lli.1",lli.1,"_lli2",lli.2, sep = "")
print(nm)

for (i in 1:nSim){
  print(i)
  temp.HH.netw<-HH.networks[[sample(1:length(HH.networks),1)]]
  epi.outbreak[[i]]<-sim.multipathogen(HH.network = temp.HH.netw, t2=t2, lambda.g = lambda.g, sigma12 = sigma12, sigma21 = sigma21, prop.immune = prop.immune, nSeeds.1 = nSeeds.1, nSeeds.2 = nSeeds.2, rho.1 = rho.1, rho.2 = rho.2, inf.path.1.h = inf.path.1.h,inf.path.1.g = inf.path.1.g, inf.path.2.h = inf.path.2.h,inf.path.2.g = inf.path.2.g, alpha.as.1=alpha.as.1,alpha.as.2=alpha.as.2, lli.1=lli.1,lli.2=lli.2, pathogen.1=pathogen.1, pathogen.2=pathogen.2, contact.reduction=contact.reduction, t.stop=t.stop, t.seed=t.seed, bc.1=bc.1, bc.2=bc.2, reinf=reinf, typeIC=typeIC)
}
scen<-paste(netw,"_nVertex",n.vertex,"_nNetw",n.networks,pathogen.1,"_&_",pathogen.2,sep ="")

name<-paste("MP_",scen,"_R1",R.1,"_R2",R.2,"_qhqg",ratio.qhqg, "_t2",t2,"_sigma12_",sigma12,"_sigma21_",sigma21,"_alpha1",alpha.as.1,"_alpha2",alpha.as.2,"_rho1",rho.1,"_rho2",rho.2,"_lli1",lli.1,"_lli2",lli.2,"_Net",netw,"_CtcRed",contact.reduction,"_PImm",prop.immune,"_tSeed",t.seed, "_bc1",bc.1,"_bc2",bc.2,".RData", sep = "")
setwd(out)
save(epi.outbreak, file = name)

