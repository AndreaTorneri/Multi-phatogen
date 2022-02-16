# 
##########################################################################
# Set of function for the co-infection script
##########################################################################

incubation.period<-function(pathogen){
  if (pathogen=="COVID-19"){
    return(rlnorm(1,meanlog = log(5.2), sdlog = log(1.7)))    
  }
  if (pathogen=="FLU-A"){
    return(rlnorm(1,meanlog = log(1.4), sdlog = log(1.51))) 
  }
#  if (pathogen=="FLU-B"){
#    return(rlnorm(1,meanlog = log(0.6), sdlog = log(1.51))) 
#  }
  if (pathogen=="RSV"){
    return(rlnorm(1,meanlog = log(4.4), sdlog = log(1.24))) 
  }
}

infectious.period.length<-function(pathogen){
  if (pathogen=="COVID-19"){
    return(15)
  }
  if (pathogen=="FLU-A"){
    return(6.25)  
  }
  if (pathogen=="FLU-B"){
    return(4.8)  
  }
  if (pathogen=="RSV"){
    return(12)  
  }
  
}

InfMeasure<-function(t,pathogen){
  if (pathogen=="COVID-19"){
    return(dgamma(t,shape = 12, rate = 2.08)/ (pgamma(15,shape = 2,rate = 2.08)*0.98))
  }
  if (pathogen=="FLU-A"){
    return(dgamma(t,shape = 3.5, rate = 1.15)/ (pgamma(6.24,shape = 3.5,rate = 1.15)*0.92))
  }
  if (pathogen=="FLU-B"){
    return(dgamma(t,shape = 3.5, rate = 1.15)/ (pgamma(6.24,shape = 3.5,rate = 1.15)*0.83))
  }
  if (pathogen=="RSV"){
    return(dgamma(t,shape = 15, rate = 2.6)/ (pgamma(12,shape = 15,rate = 2.6)) )
  }
}

long.inter.term.1<-function(t,inf.type,lli){
  #return(k.asint*(1-exp(-0.02*t)))
  return(lli)
}

long.inter.term.2<-function(t,inf.type,lli){
  #return(k.asint*(1-exp(-0.02*t)))
  return(lli)
}


#INPUT PARAMETERS:
# mu.1 (mu.2) - length infectious period for pathogen 1 (2) (remark: the infectious period is assumed to be constant over time)
# inf.path.hh.1 (inf.path.g.1) - infectiousness curve for pathogen 1 for contacts within household members (with the rest of the population). Similarly for pathogen 2
# HH.network - object of class network identifying the connection between household members. There can be only undirected links between members of the same
#             household. 
# t_2 - time at which the pathogen 2 is introduced in the population
# sigma_12 - susceptibility of acquiring 1 when having 2
# sigma_21 - susceptibility of acquiring 1 when having 2
# lambda.h (lambda.g) - rate at which individual makes contacts in a household (with member of the population)
# bc.prob.1 (bc.prob.2) - probability that an individual changes the behavior depending on the symptoms score for pathogen 1 (pathogen 2)

#OUTPUT PARAMETERS
# time events - three-column matrix identifying one of the following events(second column):
#               1.1 and infection with pathogen 1 occurs
#               1.2 and infection with pathogen 2 occurs
#               0.1 an individual recovers from an infection with pathogen 1
#               0.2 an individual recovers from an infection with pathogen 2
#               the time at which the event occur (first column) and the individual (third column)
# status matrix.1 - Matrix that represents in each row and individual, and its infectious state for pathogen 1. More precisely,
#               the first colum indicates if the individual is infected (1), susceptible (0) or recovered (NA).
#               When one individual is infected, the infection time is reported (second column) as well as the infector (third column) 
#

sim.multipathogen<-function(HH.network, t2, lambda.g, sigma21, sigma12, prop.immune, nSeeds.1,nSeeds.2, rho.1,rho.2,inf.path.1.h,inf.path.1.g,inf.path.2.h,inf.path.2.g, alpha.as.1,alpha.as.2,lli.1,lli.2, pathogen.1,pathogen.2, contact.reduction){
  
  
  
  n<-network.size(HH.network)
  hh.id<- HH.network %v% "hh_id"
  
  status.matrix.1<- data.frame(infected          = rep(0,n), # 1
                              time.of.infection = NA,        # 2
                              infector          = NA,        # 3
                              severity          = 0,         # 4 1 Symptomatic 2 Asymptomatic
                              TimeSymptomOnset  = Inf,       # 5
                              Recovery          = Inf)
  
  status.matrix.2 <-status.matrix.1
  
  
  recovery.vector.1<-rep(Inf,n) #vector giving the recovery times
  recovery.vector.2<-rep(Inf,n) #vector giving the recovery times
  
  events<-data.frame("NextCtc"=Inf, "HomeQuarantine"=Inf, "Recovery"=Inf, "NewPathogen"=Inf)
  events$NewPathogen<-t2
  
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  current.time<-0
  index.contact.within<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  index.contact.between<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  time.events<-matrix(NA,1,3)
  
  #transmission parameter dataframe: each line is an individual, the first colum is the transmsission coeffficient and the third the length of IP (needed to re-scale Viral load curve)
  transmission.parameters<-data.frame("id"=1:n,"q1h"=rep(NA,n),"q1g"=rep(NA,n),"q2h"=rep(NA,n),"q2g"=rep(NA,n), "contact_rate_within"=rep(NA,n),"contact_rate_between"=lambda.g, "susceptibility"=rep(1,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  
  for (j in 1:n){
    transmission.parameters$contact_rate_within[j]<-length(get.neighborhood(HH.network,j))
  }
  
  #Proportion of immune
  if (prop.immune>0){
    status.matrix[sample(1:n,round(prop.immune*n)),1]<--2
  }
  
  homequarantine.day.1<-rep(Inf,n)
  homequarantine.day.2<-rep(Inf,n)
  homequarantine<-rep(0,n)
  stop.quarantine<-rep(Inf,n)
  
  contact.time.within<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  contact.time.between<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  
  # first infected: randomly chosen in the population (among the susceptibles)
  first.cases<-sample(which(status.matrix.1[,1]==0),nSeeds.1)
  
  for (j in first.cases){
    first<-j
    status.matrix.1$infected[first] <- 1 
    status.matrix.1$time.of.infection[first] <- 0
    status.matrix.1$Recovery[first]<-current.time+infectious.period.length(pathogen = pathogen.1)
    if (runif(1)<rho.1){ #if symptomatic
      transmission.parameters$q1h[first]<-inf.path.1.h #A single q parameter for everyone
      transmission.parameters$q1g[first]<-inf.path.1.g #A single q parameter for everyone
      status.matrix.1$TimeSymptomOnset[first]<-current.time+incubation.period(pathogen=pathogen.1)
      homequarantine.day.1[first]<-status.matrix.1$TimeSymptomOnset[first]
      status.matrix.1$severity[first]<-1
      time.events<-rbind(time.events,c(current.time,1.1,first))
    }else{
      transmission.parameters$q1h[first]<-inf.path.1.h*alpha.as.1 #A single q parameter for everyone
      transmission.parameters$q1g[first]<-inf.path.1.g*alpha.as.1 #A single q parameter for everyone
      status.matrix.1$severity[first]<-2
      time.events<-rbind(time.events,c(current.time,1.2,first))
    }
    infectives[first]<-1
    contact.time.within$pr.ctc[first]<-ifelse(transmission.parameters$contact_rate_within[first]!=0,rexp(1,transmission.parameters$contact_rate_within[first])+current.time+status.matrix.1[first,5],Inf)       # I generate the next interarrival time for individual i
    contact.time.between$pr.ctc[first]<-rexp(1,transmission.parameters$contact_rate_between[first])+current.time+status.matrix.1[first,5] # I generate the next interarrival time for individual i
  }
  
  proposed.individual<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  recovered<-0
  err<-0
  
  #When only the first pathogen is present
  while((sum(infectives))>0 | current.time<t2){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contac
    for (i in which(index.contact.within==1) ){ # for all the individuals that has to propose a global contact
      contact.time.within$pr.ctc[i]<-rexp(1,transmission.parameters$contact_rate_within[i])+current.time# I generate the next interarrival time for individual i
      index.contact.within[i]<-0
    }
    for (i in which(index.contact.between==1) ){ # for all the individuals that has to propose a global contact
      contact.time.between$pr.ctc[i]<-rexp(1,transmission.parameters$contact_rate_between[i])+current.time# I generate the next interarrival time for individual i
      index.contact.between[i]<-0
    }
    
    contact.time.overall<-c(contact.time.within$pr.ctc, contact.time.between$pr.ctc) #overall contact times    
    recovery.vector.overall<-c(status.matrix.1$Recovery ,status.matrix.2$Recovery)
    homequarantine.day.overall<-c(homequarantine.day.1,homequarantine.day.2)
    #Phase 2: identify the next event: possible infection, recovery or the start of the new pathogen infection
    ifelse(length(which(is.na(contact.time.overall)==FALSE))>0,events$NextCtc<-min(contact.time.overall, na.rm = T),events$NextCtc<-Inf) # among all the proposed social contact between houeholds we select the minimum
    ifelse(length(which(!is.infinite(homequarantine.day.overall)))>0,events$HomeQuarantine<-min(homequarantine.day.overall),events$HomeQuarantine<-Inf ) #minimum quarantine pathogen 1
    ifelse(length(which(is.na(recovery.vector.overall)==FALSE))>0,events$Recovery<-min(recovery.vector.overall, na.rm = T),events$Recovery<-Inf) # among all the proposed social contact between houeholds we select the minimum
    
    next.evts<-colnames(events)[which(min(events)==events)]
    if (length(next.evts)>1){
      next.evts<-sample(colnames(events)[which(min(events)==events)],1)
    }
    
    if (next.evts=="NextCtc"){
      current.time<-events$NextCtc
      if (length(min(contact.time.overall, na.rm = T))>1){ #when two contacts happen at the same time
        selected.ctc<-sample(which(contact.time.overall==current.time),1) 
        if (selected.ctc!=n & selected.ctc!=2*n){
          infector<- selected.ctc %% n
        }else{
          infector<- n
        }
        if (selected.ctc<=n){
          infectee.pool<-get.neighborhood(HH.network,infector)
          if (length(infectee.pool)>1){
            infectee<-sample(infectee.pool,1) #pick a random individual in the class
          }
          if(length(infectee.pool)==1){
            infectee<-infectee.pool #pick a random individual in the class
          }
          if(length(infectee.pool)==0){
            infectee<-infector #just a trick to have acceptance rate 0 (infector is not susceptible)
          }
          index.contact.within[infector]<-1
          contact.time.within$pr.ctc[infector]<-NA
          ctc<-"hh"
        }else{
          infector<-which(contact.time.between$pr.ctc == events$NextCtc) 
          hh.members.temp<-which(hh.id==hh.id[infector])
          infectee.pool<- setdiff(1:n,hh.members.temp) #individuals not in the same class
          infectee<-sample(infectee.pool,1) #pick a random individual not in the class (and not a teacher)
          index.contact.between[infector]<-1
          contact.time.between$pr.ctc[infector]<-NA
          ctc<-"g"
        }
      }else{
        if (length(which(events$NextCtc==contact.time.within$pr.ctc))>0){ #if it is a within contact
          infector<-which(contact.time.within$pr.ctc ==events$NextCtc)
          infectee.pool<-get.neighborhood(HH.network,infector)
          if (length(infectee.pool)>1){
            infectee<-sample(infectee.pool,1) #pick a random individual in the class
          }
          if(length(infectee.pool)==1){
            infectee<-infectee.pool #pick a random individual in the class
          }
          if(length(infectee.pool)==0){
            infectee<-infector #just a trick to have acceptance rate 0 (infector is not susceptible)
          }
          index.contact.within[infector]<-1
          contact.time.within$pr.ctc[infector]<-NA
          ctc<-"hh"
        }else{
          infector<-which(contact.time.between$pr.ctc == events$NextCtc) 
          hh.members.temp<-which(hh.id==hh.id[infector])
          infectee.pool<- setdiff(1:n,hh.members.temp) #individuals not in the same class
          infectee<-sample(infectee.pool,1) #pick a random individual not in the class (and not a teacher)
          index.contact.between[infector]<-1
          contact.time.between$pr.ctc[infector]<-NA
          ctc<-"g"
        }
      }
      # compute the long and short interaction terms for pathogen.1
      if (status.matrix.2[infectee,1]==1){
        short.inter<-sigma12
      }else{
        short.inter<-1
      }
      if (status.matrix.2[infectee,1]<0){
        long.inter<-long.inter.term.2(t=status.matrix.2$Recovery[infectee],inf.type=status.matrix.2$infected[infectee], lli = lli.2)
      }else{
        long.inter<-1
      }
      
      ifelse(ctc=="g",q<-transmission.parameters$q1g[infector],q<-transmission.parameters$q1h[infector])
      acc.rate.1<-InfMeasure(t= current.time- status.matrix.1$time.of.infection[infector] ,pathogen = pathogen.1)*short.inter*long.inter*q
      if ((ctc=="g" & homequarantine[infectee]==1) | status.matrix.1$infected[infector]!=1){acc.rate.1<-0}
      if (acc.rate.1>1){err<-err+1}
      if (status.matrix.1$infected[infectee]==0 & runif(1)<acc.rate.1){
        status.matrix.1$infected[infectee] <- 1 
        status.matrix.1$time.of.infection[infectee] <- current.time
        status.matrix.1$infector[infectee] <- infector
        status.matrix.1$Recovery[infectee]<-current.time+infectious.period.length(pathogen=pathogen.1)
        if (runif(1)<rho.1){ #if symptomatic
          transmission.parameters$q1h[infectee]<-inf.path.1.h #A single q parameter for everyone
          transmission.parameters$q1g[infectee]<-inf.path.1.g #A single q parameter for everyone
          status.matrix.1$TimeSymptomOnset[infectee]<-current.time+incubation.period(pathogen=pathogen.1)
          homequarantine.day.1[infectee]<-status.matrix.1$TimeSymptomOnset[infectee]
          status.matrix.1$severity[infectee]<-1
          time.events<-rbind(time.events,c(current.time,1.1,infectee))
        }else{
          transmission.parameters$q1g[infectee]<-inf.path.1.g*alpha.as.1 #A single q parameter for everyone
          transmission.parameters$q1h[infectee]<-inf.path.1.h*alpha.as.1 #A single q parameter for everyone
          status.matrix.1$severity[infectee]<-2
          time.events<-rbind(time.events,c(current.time,1.0,infectee))
        }
        if (infectives[infectee]==0){
          infectives[infectee]<-1
          contact.time.within$pr.ctc[infectee]<-ifelse(transmission.parameters$contact_rate_within[infectee]!=0,rexp(1,transmission.parameters$contact_rate_within[infectee])+current.time,Inf)       # I generate the next interarrival time for individual i
          if (homequarantine[infectee]==0){
            contact.time.between$pr.ctc[infectee]<-rexp(1,transmission.parameters$contact_rate_between[infectee])+current.time # I generate the next interarrival time for individual i
          }
        }
      }
      
      
      # compute the long and short interaction terms for pathogen.2
      
      if (status.matrix.1$infected[infectee]==1){
        short.inter<-sigma21
      }else{
        short.inter<-1
      }
      if (status.matrix.1$infected[infectee]<0){
        long.inter<-long.inter.term.1(t=status.matrix.1$Recovery[infectee],inf.type=status.matrix.1$infected[infectee],lli = lli.1)
      }else{
        long.inter<-1
      }
      ifelse(ctc=="g",q<-transmission.parameters$q2g[infector],q<-transmission.parameters$q2h[infector])
      acc.rate.2<-InfMeasure(t=(current.time-status.matrix.2$time.of.infection[infector]),pathogen = pathogen.2)*short.inter*long.inter*q
      if ((ctc=="g" & homequarantine[infectee]==1) | status.matrix.2$infected[infector]!=1){acc.rate.2<-0}
      if (acc.rate.2>1){err<-err+1}
      if (status.matrix.2$infected[infectee]==0 & runif(1)<acc.rate.2){
        status.matrix.2$infected[infectee] <- 1 
        status.matrix.2$time.of.infection[infectee] <- current.time
        status.matrix.2$infector[infectee] <- infector
        status.matrix.2$Recovery[infectee]<-current.time+infectious.period.length(pathogen=pathogen.2)
        if (runif(1)<rho.2){ #if symptomatic
          transmission.parameters$q2h[infectee]<-inf.path.2.h #A single q parameter for everyone
          transmission.parameters$q2g[infectee]<-inf.path.2.g #A single q parameter for everyone
          status.matrix.2$TimeSymptomOnset[infectee]<-current.time+incubation.period(pathogen=pathogen.2)
          homequarantine.day.2[infectee]<-status.matrix.2$TimeSymptomOnset[infectee]
          status.matrix.2$severity[infectee]<-1
          time.events<-rbind(time.events,c(current.time,2.1,infectee))
        }else{
          transmission.parameters$q2g[infectee]<-inf.path.2.g*alpha.as.2 #A single q parameter for everyone
          transmission.parameters$q2h[infectee]<-inf.path.2.h*alpha.as.2 #A single q parameter for everyone
          status.matrix.2$severity[infectee]<-2
          time.events<-rbind(time.events,c(current.time,2.2,infectee))
        }
        if (infectives[infectee]==0){
          infectives[infectee]<-1
          contact.time.within$pr.ctc[infectee]<-ifelse(transmission.parameters$contact_rate_within[infectee]!=0,rexp(1,transmission.parameters$contact_rate_within[infectee])+current.time,Inf)       # I generate the next interarrival time for individual i
          if (homequarantine[infectee]==0){
            contact.time.between$pr.ctc[infectee]<-rexp(1,transmission.parameters$contact_rate_between[infectee])+current.time # I generate the next interarrival time for individual i
          }
        }
      }
    }
    
    
    if (next.evts=="HomeQuarantine"){
      current.time<-events$HomeQuarantine
      quarantined.individuals<-which(homequarantine.day.overall==current.time)
      for (k in quarantined.individuals){
        if (k != n & k!= 2*n) {
          temp.ind<- k %% n          
        }else{
          temp.ind<-n
        }
        if (k>n){
          homequarantine.day.2[temp.ind]<-Inf
          stop.quarantine[temp.ind]<-status.matrix.2$Recovery[temp.ind]
        }else{
          homequarantine.day.1[temp.ind]<-Inf
          stop.quarantine[temp.ind]<-status.matrix.1$Recovery[temp.ind]
        }
        if (homequarantine[temp.ind]==1){ #individual is already in quarantine for the other disease
          stop.quarantine[temp.ind]<-max(status.matrix.2$Recovery[temp.ind],status.matrix.1$Recovery[temp.ind])
        }
        homequarantine[temp.ind]<-1
        contact.time.between$pr.ctc[temp.ind]<-NA
        index.contact.between[temp.ind]<-0
        transmission.parameters$contact_rate_within[temp.ind]<-transmission.parameters$contact_rate_within[temp.ind]*contact.reduction
      }
    }
    
    
    if (next.evts=="Recovery"){
      current.time<-events$Recovery
      temp.recovered<-which(recovery.vector.overall==events$Recovery)
      for (recovered in temp.recovered){
        if (recovered!= n & recovered!=n*2){
          if (recovered > n){
            recovered<- recovered %% n        
            status.matrix.2$infected[recovered]<--1
            status.matrix.2$Recovery[recovered]<-Inf
            time.events<-rbind(time.events,c(current.time,-2,recovered))
            if (status.matrix.1$infected[recovered]!=1){
              infectives[recovered]<-0
              contact.time.between$pr.ctc[recovered]<-NA
              contact.time.within$pr.ctc[recovered]<-NA
              index.contact.within[recovered]<-0
              index.contact.between[recovered]<-0
            }else{
              if (homequarantine[recovered]==0){
                index.contact.between[recovered]<-1
              }
            }
          }else{
            status.matrix.1$infected[recovered]<--1
            status.matrix.1$Recovery[recovered]<-Inf
            time.events<-rbind(time.events,c(current.time,-1,recovered))
            if (status.matrix.2$infected[recovered]!=1){
              infectives[recovered]<-0
              contact.time.between$pr.ctc[recovered]<-NA
              contact.time.within$pr.ctc[recovered]<-NA
              index.contact.within[recovered]<-0
              index.contact.between[recovered]<-0
            }else{
              if (homequarantine[recovered]==0){
                index.contact.between[recovered]<-1
              }
            }
          }
        }else{
          if (recovered == 2*n){
            recovered<- n        
            status.matrix.2$infected[recovered]<--1
            status.matrix.2$Recovery[recovered]<-Inf
            time.events<-rbind(time.events,c(current.time,-2,recovered))
            if (status.matrix.1$infected[recovered]!=1){
              infectives[recovered]<-0
              contact.time.between$pr.ctc[recovered]<-NA
              contact.time.within$pr.ctc[recovered]<-NA
              index.contact.within[recovered]<-0
              index.contact.between[recovered]<-0
            }else{
              if (homequarantine[recovered]==0){
                index.contact.between[recovered]<-1
              }
            }
          }else{
            status.matrix.1$infected[recovered]<--1
            status.matrix.1$Recovery[recovered]<-Inf
            time.events<-rbind(time.events,c(current.time,-1,recovered))
            if (status.matrix.2$infected[recovered]!=1){
              infectives[recovered]<-0
              contact.time.between$pr.ctc[recovered]<-NA
              contact.time.within$pr.ctc[recovered]<-NA
              index.contact.within[recovered]<-0
              index.contact.between[recovered]<-0
            }else{
              if (homequarantine[recovered]==0){
                index.contact.between[recovered]<-1
              }
            }
          }
        }
        if (stop.quarantine[recovered]==current.time){
          homequarantine[recovered]<-0
          stop.quarantine[recovered]<-Inf
        }
      }
    }
    
    
    if (next.evts=="NewPathogen"){
      current.time<-events$NewPathogen
      events$NewPathogen<-Inf
      first.cases<-sample(1:n,nSeeds.2)
      for (j in first.cases){
        first<-j
        status.matrix.2$infected[first] <- 1 
        status.matrix.2$time.of.infection[first] <- current.time
        status.matrix.2$Recovery[first]<-current.time+infectious.period.length(pathogen=pathogen.2)
        if (runif(1)<rho.2){ #if symptomatic
          transmission.parameters$q2h[first]<-inf.path.2.h #A single q parameter for everyone
          transmission.parameters$q2g[first]<-inf.path.2.g #A single q parameter for everyone
          status.matrix.2$severity[first]<-1
          status.matrix.2$TimeSymptomOnset[first]<-current.time+incubation.period(pathogen=pathogen.2)
          homequarantine.day.2[first]<-status.matrix.2$TimeSymptomOnset[first]
          time.events<-rbind(time.events,c(current.time,2.1,first))
        }else{
          transmission.parameters$q2h[first]<-inf.path.2.h*alpha.as.1 #A single q parameter for everyone
          transmission.parameters$q2g[first]<-inf.path.2.g*alpha.as.1 #A single q parameter for everyone
          status.matrix.2$severity[first]<-2
          time.events<-rbind(time.events,c(current.time,2.0,first))
        }
        if (infectives[first]==0){
          infectives[first]<-1
          contact.time.within$pr.ctc[first]<-ifelse(transmission.parameters$contact_rate_within[first]!=0,rexp(1,transmission.parameters$contact_rate_within[first])+current.time,Inf)       # I generate the next interarrival time for individual i
          if (homequarantine[first]==0){
            contact.time.between$pr.ctc[first]<-rexp(1,transmission.parameters$contact_rate_between[first])+current.time # I generate the next interarrival time for individual i
          }
        }
      }
    }
    
  }
  
  #When also the other pathogen is present.
  time.events<-time.events[-1,]
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  first.cases.1<-which(status.matrix.1$time.of.infection==0)
  for (o in first.cases.1){
    temp.sec.cases<-NULL
    ifelse(length(which(status.matrix.1$infector ==o)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix.1$infector==o))),temp.sec.cases<-c(temp.sec.cases,0))
  }
  Rt1<-mean(temp.sec.cases)
  
  first.cases.2<-which(status.matrix.2$time.of.infection==0)
  for (o in first.cases.2){
    temp.sec.cases<-NULL
    ifelse(length(which(status.matrix.2$infector==o)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix.2$infector==o))),temp.sec.cases<-c(temp.sec.cases,0))
  }
  Rt2<-mean(temp.sec.cases)
  C1<-nSeeds.1
  C2<-nSeeds.2
  Y1<-nSeeds.1
  Y2<-nSeeds.2
  last.day<-round(max(time.events[,1]))
  
  for (i in 1:last.day){
    temp.time<-setdiff(which(time.events[,1]>i),which(time.events[,1]>i+1))
    temp.inf.1<-c(which(time.events[temp.time,2]==1.0),which(time.events[temp.time,2]==1.1))
    temp.inf.2<-c(which(time.events[temp.time,2]==2.0),which(time.events[temp.time,2]==2.1))
    
    temp.time.1<-setdiff(1:length(time.events[,1]),which(time.events[,1]>i+1))
    C1<- c(C1,length((which(time.events[temp.time.1,2]==1.0)))+length((which(time.events[temp.time.1,2]==1.1)))-length((which(time.events[temp.time.1,2]==-1))))
    C2<- c(C2,length((which(time.events[temp.time.1,2]==2.0)))+length((which(time.events[temp.time.1,2]==2.1)))-length((which(time.events[temp.time.1,2]==-2))))
    Y1<- c(Y1,length(temp.inf.1))
    Y2<-c(Y2,length(temp.inf.2))
    
    if (length(temp.inf.1)>0){
      newly.infected<-time.events[temp.time[temp.inf.1],3]
      temp.sec.cases<-NULL
      for (k in newly.infected) {
        ifelse(length(which(status.matrix.1$infector==k)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix.1$infector==k))),temp.sec.cases<-c(temp.sec.cases,0))
      }
      Rt1<-c(Rt1,mean(temp.sec.cases))
    }else{
      Rt1<-c(Rt1,NA)
    }
    if (length(temp.inf.2)>0){
      newly.infected<-time.events[temp.time[temp.inf.2],3]
      temp.sec.cases<-NULL
      for (k in newly.infected) {
        ifelse(length(which(status.matrix.2$infector==k)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix.2$infector==k))),temp.sec.cases<-c(temp.sec.cases,0))
      }
      Rt2<-c(Rt2,mean(temp.sec.cases))
    }else{
      Rt2<-c(Rt2,NA)
    }
  }
  epi.details<-data.frame("Days"=0:last.day, "Incidence1"=Y1,"Incidence2"=Y2, "Prevalence1"=C1,"Prevalence2"=C2,"Rt1"=Rt1,"Rt2"=Rt2)
  FinalSize<-data.frame("FinalSize1"=length(which(status.matrix.1$infected==-1)),"FinalSize2"=length(which(status.matrix.2$infected==-1)))
  PeakIncidence<-data.frame("PeakIncidence1"=max(epi.details$Incidence1),"TimePeakIncidence1"=which(epi.details$Incidence1==max(epi.details$Incidence1))[1],"PeakIncidence2"=max(epi.details$Incidence2),"TimePeakIncidence2"=which(epi.details$Incidence2==max(epi.details$Incidence2))[1] )
  PeakPrevalence<-data.frame("PeakPrevalence1"=max(epi.details$Prevalence1),"TimePeakPrevalence1"=which(epi.details$Prevalence1==max(epi.details$Prevalence))[1],"PeakPrevalence2"=max(epi.details$Prevalence2),"TimePeakPrevalence2"=which(epi.details$Prevalence2==max(epi.details$Prevalence2))[1] )
  return(list(time.events=time.events, status.matrix.1=status.matrix.1, status.matrix.2=status.matrix.2,epi.details=epi.details, FinalSize=FinalSize, PeakIncidence=PeakIncidence, PeakPrevalence=PeakPrevalence))
}














