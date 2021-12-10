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
R.1= as.numeric(args[13]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R.1=",R.1)
R.2= as.numeric(args[14]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",R.2=",R.2)
lli.k= as.numeric(args[15]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",lli.k=",lli.k)
pathogen.1= as.character(args[16]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",pathogen.1=",pathogen.1)
pathogen.2= as.character(args[17]) # short-term interaction parameter: acquiring 1 while having 2 (this value can)
cat(",pathogen.2=",pathogen.2)


#Input parameters - fixed
library(ergm)
library(RGeode)
#load("HH_complete.RData")
load("sim_basis_complete_n_1000.RData")

HH.networks<-HH_sim
lambda.h<-2.9 #(computed as the average number of links individuals have within households)
Net<-"HH_Eva_ERGM_qh8.27qg"

#ERGM.EGO
#Networks<-readRDS("~/Documents/Work/PhD/Co-infection/Script/sim_EgoHHcombined_basis_complete_n_10.rds")
#HH.network<-Networks[[1]]
#id<- HH.network %v% ".NetworkID"
#set.vertex.attribute(HH.network, "hh_id", value = id)
#n.link<-NULL
#hh_size<-NULL
#for (i in 1:length(id)) {
#  n.link<-c(n.link,length(get.neighborhood(HH.network,i,type = "out")))
#  hh_size<-c(hh_size, which(id==i))
#}
#set.vertex.attribute(HH.network, "hh_size", value = hh_size)
#avg.ct<-mean(n.link)
#lambda.h<-avg.ct
lambda.g<-8.29

#SettingInfectivityParameters
# direct computation
source("R_comp_netw.R")
ratio_hhgl<-8.27*lambda.h/lambda.g
R.rif<-R.1
nSim<-100
tol<-0.05
nSeed<-3082021
set.seed(nSeed)
trs.prms<-R0.comp(ratio_hhgl=ratio_hhgl, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif,prob.asym=(1-rho.1),asymp.rel.inf=alpha.as.1,lambda.h = lambda.h)

#load data
inf.path.1.h<-trs.prms$beta.h/lambda.h
inf.path.1.g<-trs.prms$beta.g/lambda.g

R.rif<-R.2
trs.prms<-R0.comp(ratio_hhgl=ratio_hhgl, HH.network = HH.networks, nSim = nSim, tol=tol,R.rif = R.rif, prob.asym=(1-rho.2),asymp.rel.inf=alpha.as.2,lambda.h = lambda.h)

inf.path.2.h<-trs.prms$beta.h/lambda.h
inf.path.2.g<-trs.prms$beta.g/lambda.g

#Compute the reproduction number related to the selected network. 
source("function.multipathogen.new.R")
nSim<-2000
epi.outbreak<-list()
nSeed<-1062021
set.seed(nSeed)

nm<-paste("t2_",t2, "_sigma12_",sigma12,"_sigma21_",sigma21,"_qh1_",inf.path.1.h,"_qg1_",inf.path.1.g,"_qh2_",inf.path.2.h,"_qg2_",inf.path.2.g, "_rho1_",rho.1,"_rho2_",rho.2,"_alpha1_",alpha.as.1,"_alpha2_",alpha.as.2,"_Path1",pathogen.1,"_Path2",pathogen.2,"_lli.k",lli.k, sep = "")
print(nm)

for (i in 1:nSim){
  print(i)
  temp.HH.netw<-HH.networks[[sample(1:length(HH.networks),1)]]
  epi.outbreak[[i]]<-sim.multipathogen(HH.network = temp.HH.netw, t2=t2, lambda.g = lambda.g, sigma12 = sigma12, sigma21 = sigma21, prop.immune = prop.immune, nSeeds.1 = nSeeds.1, nSeeds.2 = nSeeds.2, rho.1 = rho.1, rho.2 = rho.2, inf.path.1.h = inf.path.1.h,inf.path.1.g = inf.path.1.g, inf.path.2.h = inf.path.2.h,inf.path.2.g = inf.path.2.g, alpha.as.1=alpha.as.1,alpha.as.2=alpha.as.2, lli.k=lli.k, pathogen.1=pathogen.1, pathogen.2=pathogen.2)
  }

scen<-paste(pathogen.1,"_&_",pathogen.2,sep ="")

name<-paste("MP_",scen,"_R1",R.1,"_R2",R.2,"_t2",t2,"_sigma12_",sigma12,"_sigma21_",sigma21,"_alpha1",alpha.as.1,"_alpha2",alpha.as.2,"_rho1",rho.1,"_rho2",rho.2,"_lli",lli.k,"_Net",Net,"_.RData", sep = "")
setwd(out)
save(epi.outbreak, file = name)

