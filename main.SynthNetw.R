args <- commandArgs(trailingOnly = TRUE)
out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
n.vertex = as.numeric(args[3]) #time at which pathogen 2 is insert in the population
cat(",n.vertex=",n.vertex)
n.networks = as.numeric(args[4]) #time at which pathogen 2 is insert in the population
cat(",n.networks=",n.networks)

source("CreatingSynthNetw.R")
load("DensityHHsize.RData")

HH.networks<-creat.synth.netw(n.vertex = n.vertex,n.networks = n.networks,density_by_hh_size = density_by_hh_size)
name<-paste("HH_Networks","_nVertex",n.vertex,"_nNetw",n.networks,".RData",sep = "")
save(HH.networks,file = name)