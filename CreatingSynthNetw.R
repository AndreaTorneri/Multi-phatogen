
##################################################
### Function to create a new synthetic network ###
### Called by the script 'main.SynthNetw.R'    ###
##################################################

## This function generates a given number (n.networks) of new synthetic networks 
## with a given number of households (n.vertex) based on a given distribution of 
## household sizes (density_by_hh_size). 

create.synth.netw <- function(n.vertex, n.networks, density_by_hh_size){
  # load the needed library
  library(network)
  
  # Create an empty list to store the networks
  HH.netw <- list()
  # determine a maximal household size according to this maximal household size
  max.hhsz <- 7 
  # extract the relevant rows from the household density data set and compute the relative proportions
  data.hhnetw <- density_by_hh_size[1:max.hhsz,]
  prob.hhsize <- data.hhnetw$household_size_freq/sum(data.hhnetw$household_size_freq)
  
  # this count is needed to keep track of how large the network is during the creation process
  count <- 0
  
  # create a data frame to store the frequency of households of each size
  # the first column contains the household sizes (from 1 to max.hhsz) and the second column the frequency
  hh.size.synth <- data.frame("Size" = 1:max.hhsz, "Freq" = 0)
  
  ## sample household sizes between 1 and max.hhsz with according probabilities and keep track of the numbers 
  ## of households of specific sizes in hh.size.synth
  for (w in 1:n.networks){ # for each network to create
    # create households of different sizes until the network is of desired size. 
    while(count<n.vertex){ # check whether the total population size is not yet exceeded
      # k is the household size (sampled)
      k <- sample(1:max.hhsz,1,prob = prob.hhsize)
      # update the data.frame that keeps track of frequency of networks of different sizes
      hh.size.synth$Freq[which(hh.size.synth$Size==k)] <- hh.size.synth$Freq[which(hh.size.synth$Size==k)]+1
      # update the count
      count <- count+k
    }
    
    # create a square matrix of size nxn to be the adjacency matrix of network members
    adj.mat <- matrix(data=0,nrow = count, ncol = count)
    # create a vector of length n
    hh.size.vect <- rep(1,count)
    
    ## fill in hh.size.vect
    temp2 <- hh.size.synth$Freq[1] * (1)
    for(i in 2:max.hhsz){ # for all household sizes starting from size 2
      if(hh.size.synth[i,2]!=0){ # extra check for the case that a certain household size is not represented
        temp1 = temp2
        temp2 = hh.size.synth$Freq[i] * (i) + temp1
        hh.size.vect[(temp1+1):(temp2)] = i
        # update the adjacency matrix
        for(j in seq(temp1+1, temp2, by = i)){
          household = seq(j,j+i-1,1)
          combinations = combn(household,2)
          for(x in 1:ncol(combinations)){
            adj.mat[combinations[1,x], combinations[2,x]] = i
            adj.mat[combinations[2,x], combinations[1,x]] = i
          }
        }
      }
    }
    
    # Constructing hh.id, a vector that contains the household size of each individual at its index. 
    hh.id.vect <- c(1:hh.size.synth$Freq[1]) # fill vector for household size = 1
    for(i in 2:max.hhsz){ # idem for household sizes > 1
      if(hh.size.synth[i,2] != 0){
        temp.seq <- seq(hh.id.vect[length(hh.id.vect)]+1, hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[i]) # for all households of size 2
        hh.id.vect <- c(hh.id.vect, rep(temp.seq,each = i))
      }
    }
    
    vert.attr<-list()
    vert.attr.name<-list()
    vert.attr[[1]]<-hh.size.vect
    vert.attr[[2]]<-hh.id.vect
    vert.attr.name[[1]]<-"hh_size"
    vert.attr.name[[2]]<-"hh_id"
    HH.netw[[w]]<-network(adj.mat, 
                          directed = FALSE, 
                          vertex.attr = vert.attr, 
                          vertex.attrnames = vert.attr.name)
  }
  return(HH.netw)
}



