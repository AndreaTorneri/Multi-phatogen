###########################################################################
### This script creates a list of synthetic networks                   ####
### It relies on functions defined in the script 'CreatingSynthNetw.R' ####
###########################################################################

### These parameters come from an external file 'CreatingNetw.csv' 
args <- commandArgs(trailingOnly = TRUE)
# working directory
out <- args[1] 
cat(",out=",out)
# number of cores to run in parallel
cores <- as.numeric(args[2]) 
cat(",cores=",cores)
# number of individuals in the network
n.vertex <- as.numeric(args[3]) 
cat(",n.vertex=",n.vertex)
# numbers of networks to create
n.networks <- as.numeric(args[4]) 
cat(",n.networks=",n.networks)

### For testing purposes, create a small network with the following parameters:
n.vertex <- 100
n.networks <- 1

### Load necessary files
## This script relies on functions defined in the script 'CreatingSynthNetw.R'
source("CreatingSynthNetw.R")
## Dataset containing densities for all household sizes. 
# Set working directory if necessary
# setwd("C:/Users/LUCP13441/Documents/GitHub/Multi-phatogen")
load("DensityHHsize.RData")

# Call the function to create a new list of networks
HH.networks <- create.synth.netw(n.vertex = n.vertex,
                                 n.networks = n.networks,
                                 density_by_hh_size = density_by_hh_size)
# Determine the name of the file and save it
name <- paste("HH_Networks", "_nVertex", n.vertex, "_nNetw", n.networks, ".RData", sep = "")
save(HH.networks, file = name)
