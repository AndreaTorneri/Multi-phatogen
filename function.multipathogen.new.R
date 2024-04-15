# 
##########################################################################
# Set of function for the co-infection script
##########################################################################

incubation.period<-function(pathogen){
  if (pathogen=="COVID-19" | pathogen == "DELTA" | pathogen== "OMICRON"){
    #return(rlnorm(1,meanlog = log(5.2), sdlog = log(1.7)))
    return(5.288462)
  }
  if (pathogen=="FLU-A"){
    #return(rlnorm(1,meanlog = log(1.4), sdlog = log(1.51))) 
    return(2)
  }
  if (pathogen=="RSV"){
    return(rlnorm(1,meanlog = log(4.4), sdlog = log(1.24))) 
  }
  if (pathogen == "XP"){
    return(1)
  }
  if (pathogen == "XS"){
    return(2)
  }
  if (pathogen == "XA"){
    return(3)
  }
}

infectious.period.length<-function(pathogen){
  if (pathogen=="COVID-19" | pathogen=="OMICRON" | pathogen=="DELTA"){
    return(12)
  }
  if (pathogen=="FLU-A"){
    return(8)  
  }
  if (pathogen=="FLU-B"){
    return(4.8)  
  }
  if (pathogen == "XP" | pathogen == "XS" | pathogen == "XA"){
      return(4)
  }
}

InfMeasure<-function(t,pathogen){
  if (pathogen=="COVID-19" | pathogen=="DELTA" | pathogen=="OMICRON"){
    return(dgamma(t,shape = 12, rate = 2.08)/ (pgamma(12,shape = 12,rate = 2.08)))
  }
  if (pathogen=="FLU-A"){
    # Setting infectiousness measure according to Carrat et al. (2008) for H1N1
    #VL<-data.frame(x=0:8,y=c(0,1.75,3,2.5,1.8,1.25,0.75,0.5,0)) 
    #vl.flu<-nlsLM(y~a*dgamma(x=x,shape = s1,scale = sc1),start = list(a=10,s1=1.5,sc1=1.5),data = VL, weights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,100))
    #prms<-vl.flu$m$getPars()
    #print(prms)
    #plot(seq(0,8,0.1),prms[1]*dgamma(x=seq(0,8,0.1),shape = prms[2],scale = prms[3]),ylim = c(0,3.1))
    #points(VL, col="red")
    #f.vl<-function(t){
    #  return(prms[1]*dgamma(t,shape = prms[2],scale = prms[3]))
    #}
    #k<-integrate(f.vl,lower = 0,upper = 8)
    #prms[1]/ k$value
    return(1.001592*dgamma(t,shape = 4.604016,scale = 0.5922317))
    
    
    #return(dgamma(t,shape = 3.5, rate = 1.15)/ (pgamma(6.24,shape = 3.5,rate = 1.15)))
  }
  if (pathogen=="FLU-B"){
    return(dgamma(t,shape = 3.5, rate = 1.15)/ (pgamma(6.24,shape = 3.5,rate = 1.15)))
  }
  if (pathogen=="RSV"){
    return(dgamma(t,shape = 15, rate = 2.6)/ (pgamma(12,shape = 15,rate = 2.6)) )
  }
  if (pathogen == "XP" | pathogen == "XS" | pathogen == "XA"){
    if (t <=2){
      return(0.25*t)
    }else{
      return(-0.25*t+1)
    }
  }
}


VaccineEffectiveness<-function(t,typeIC){
  #Curve approximating symptomatic infection Omicron Qatar Chemateilly et al. 2022
  if (typeIC==1){
    return(6104.4743*dlnorm(t,meanlog = 4.3125, sdlog = 0.9887))
  }
  if (typeIC==2){
    return(4698.209*dgamma(t,shape = 2.026, scale = 32.904))
  }
  if (typeIC==3){
    return(7568.209*dlnorm(t,meanlog = 4.599, sdlog = 1.118))
  }
}

VE.flu<-function(){
  return((1-0.7))
}

VE.COVID<- function() {
  return(1-0.88)
}



# LLImmlev<-function(pathogen.v1,pathogen.v2,status.matrix.v1,status.matrix.v2,infectee,lli,current.time,reinf,typeIC,het.vac){ #pathogen.v1 is the infection the infectee might catch
#   value<-1
#   
#   if (reinf==0){
#     if (status.matrix.v1$infected[infectee]!=0){
#       value<-0
#     }
#     if (status.matrix.v1$infected[infectee]==0 & status.matrix.v2$infected[infectee]==-1){
#       value<-lli
#     }
#   }else{
#     if (pathogen.v1 =="FLU-A"){ # For FLU we assume that re-infection is possible only for a vaccine-induced immunity
#       if (status.matrix.v1$Immunity[infectee]==1){
#         if (status.matrix.v1$infected[infectee]==0){
#           value<-VE.flu()
#         }else{
#           value<-0
#         }
#       }else{
#         if (status.matrix.v1$infected[infectee]==0){
#           value<-1
#           if (status.matrix.v2$infected[infectee]==-1){
#             value<-lli
#           }
#         }else{
#           value<-0
#         }
#       }
#     }     
#     if (pathogen.v1 =="COVID-19"){
#       if (status.matrix.v1$Immunity[infectee]==1){
#         value<-0
#        if (status.matrix.v1$infected[infectee]==0 | status.matrix.v1$infected[infectee]==-1){
#           value<-VE.COVID()
#        }
#       #  if (status.matrix.v1$infected[infectee]==0 | status.matrix.v1$infected[infectee]==-1){
#       #    value<-VE.COVID()
#       #  }else{
#       #    value<-0
#       #  }
#         
#       }else{
#         value<-0
#         if (status.matrix.v1$infected[infectee]==0){
#           value<-1
#           if (status.matrix.v2$infected[infectee]==-1){
#             value<-lli
#           }
#           if (status.matrix.v2$Immunity[infectee]==1 & het.vac==1){ #heterologous effect vaccination
#             value<-(1-0.297) #Tayar et al. 2023
#           }
#         }
#         if (status.matrix.v1$infected[infectee]==-1){
#           value<-VE.COVID()
#         }
#       }
#     }     
#     
#     if (pathogen.v1 =="DELTA"){
#       if (status.matrix.v1$Immunity[infectee]==1){
#         if (status.matrix.v1$infected[infectee]==0 | status.matrix.v1$infected[infectee]==-1){
#           value<-VE.COVID()
#         }else{
#           value<-0
#         }
#       }else{
#         value<-0
#         if (status.matrix.v1$infected[infectee]==0){
#           if (status.matrix.v2$infected[infectee]==-1){
#             value<-1-(VaccineEffectiveness(t=status.matrix.v2$time.of.infection[infectee],typeIC = typeIC))/100
#           }else{
#             value<-1
#           }
#         }
#         if (status.matrix.v1$infected[infectee]==-1){
#           value<-VE.COVID()
#         }
#       }
#     }     
#     
#     if (pathogen.v1 =="OMICRON"){
#       
#       if (status.matrix.v1$Immunity[infectee]==1){
#         value<-0
#         if (status.matrix.v1$infected[infectee]==0){
#           if (status.matrix.v2$infected[infectee]==-1){
#             value<-1-(VaccineEffectiveness(t=status.matrix.v2$time.of.infection[infectee],typeIC = typeIC))/100
#           }else{
#             value<-1-(VaccineEffectiveness(t=0,typeIC = typeIC))/100
#           }
#         }
#         if (status.matrix.v1$infected[infectee]==-1){
#             value<-VE.COVID()
#         }
#       }else{
#         value<-0
#         if (status.matrix.v1$infected[infectee]==0){
#           if (status.matrix.v2$infected[infectee]==-1){
#             value<-1-(VaccineEffectiveness(t=status.matrix.v2$time.of.infection[infectee],typeIC = typeIC))/100
#           }else{
#             value<-1
#           }
#         }
#         if (status.matrix.v1$infected[infectee]==-1){
#           value<-VE.COVID()
#         }
#       }
#     }     
#     
#   }
#   return(value)
# }





LLImmlev.basic<-function(status.matrix.v2,infectee,lli,current.time,typeIC,t.imm.lim,pathogen1,pathogen2){ #pathogen.v1 is the infection the infectee might catch
  value<-1
  if (status.matrix.v2$infected[infectee]==-1){
    #t.sinc.inf<-current.time-(status.matrix.v2$Recovery[infectee]+infectious.period.length(pathogen = pathogen2))
    t.sinc.inf<-current.time-(status.matrix.v2$time.of.infection[infectee]+infectious.period.length(pathogen = pathogen2)) #here we are assuming a constant IPL
    if (t.sinc.inf<t.imm.lim){
      if (typeIC==1){
        value<-lli
      }
      if (typeIC==2){
        value<-(t.sinc.inf/t.imm.lim)
      }
      if (typeIC==3){
        value<-0.25+(t.sinc.inf/(2*t.imm.lim))
        }
      if (typeIC==4){
#        if (pathogen1=="COVID-19"){
          value<-2.5-t.sinc.inf/10
#        }else{
#          if (t.sinc.inf<10){
#            value<-(t.sinc.inf/10)
#          }
#        }
      }
      if (typeIC==5){
#        if (pathogen1=="COVID-19"){
          value<-lli
#        }else{
#          if (t.sinc.inf<10){
#           value<-(t.sinc.inf/10)
#          }
#        }
      }
      if (typeIC==6){
#        if (pathogen1=="COVID-19"){
          value<-3-t.sinc.inf/5
#       }else{
#          if (t.sinc.inf<10){
#            value<-(t.sinc.inf/10)
#          }
#        }
      }
    }  
  }
  return(value)
}


#Computing summary measures
# Rt
comp.RT<-function(status.matrix,individual,Rt){
  infectees<-which(status.matrix$infector==individual)
  Rt.temp<-0
  if (length(infectees)>0){
    for (i1 in 1:length(infectees)){
      if (status.matrix$time.of.infection[infectees[i1]]>status.matrix$time.of.infection[individual]){
        Rt.temp<-Rt.temp+1
      }
    }
  }
  Rt<-rbind(Rt,c(status.matrix$time.of.infection[individual],Rt.temp))
  return(Rt)
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

sim.multipathogen<-function(HH.network, t2, lambda.g, sigma21, sigma12, prop.immune, nSeeds.1,nSeeds.2, rho.1,rho.2,inf.path.1.h,inf.path.1.g,inf.path.2.h,inf.path.2.g, alpha.as.1,alpha.as.2,lli.1,lli.2, pathogen.1,pathogen.2, contact.reduction,t.stop,t.seed, bc.1,bc.2,reinf,typeIC, het.vac, t.imm.lim){
  
  
  
  
  n<-network.size(HH.network)
  hh.id<- HH.network %v% "hh_id"
  hh.size<- HH.network %v% "hh_size"
  
  status.matrix.1<- data.frame(infected          = rep(0,n),  # 1
                              time.of.infection  = NA,        # 2
                              infector           = NA,        # 3
                              severity           = 0,         # 4 1 Symptomatic 2 Asymptomatic
                              TimeSymptomOnset   = Inf,       # 5
                              Immunity           = 0,         # 0 no immunity, 1 vaccinated
                              Recovery           = Inf)
  
  status.matrix.2 <-status.matrix.1
  recovery.vector.1<-rep(Inf,n) #vector giving the recovery times
  recovery.vector.2<-rep(Inf,n) #vector giving the recovery times
  
  events<-data.frame("NextCtc"=Inf, "HomeQuarantine"=Inf, "Recovery"=Inf, "NewPathogen"=Inf, "NewSeeding1"=Inf, "NewSeeding2"=Inf )
  events$NewPathogen<-t2
  events$NewSeeding1<-t.seed
  events$NewSeeding2<-t2+t.seed

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
    if (pathogen.1=="DELTA" & pathogen.2=="OMICRON"){
      immuned.individuals<-sample(1:n,round(prop.immune*n))
      status.matrix.1$Immunity[immuned.individuals]<-1
      status.matrix.2$Immunity<-status.matrix.1$Immunity
    }else{
      immuned.individuals<-sample(1:n,round(prop.immune*n))
      status.matrix.1$Immunity[immuned.individuals]<-1
      immuned.individuals<-sample(1:n,round(prop.immune*n))
      status.matrix.2$Immunity[immuned.individuals]<-1
    }
  }
  
  
  homequarantine.day.1<-rep(Inf,n)
  homequarantine.day.2<-rep(Inf,n)
  homequarantine<-rep(0,n)
  stop.quarantine<-rep(Inf,n)
  
  contact.time.within<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  contact.time.between<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  
  # first infected: randomly chosen in the population (among the susceptibles)
  potential.seeds<-which(hh.size==2)
  first.cases<-sample(potential.seeds[which(status.matrix.1[potential.seeds,1]==0)],nSeeds.1)
  
  for (j in first.cases){
    first<-j
    status.matrix.1$infected[first] <- 1 
    status.matrix.1$time.of.infection[first] <- 0
    status.matrix.1$Recovery[first]<-current.time+infectious.period.length(pathogen = pathogen.1)
   # if (runif(1)<rho.1){ #if symptomatic #index cases are always symptomatic individuals
      transmission.parameters$q1h[first]<-inf.path.1.h #A single q parameter for everyone
      transmission.parameters$q1g[first]<-inf.path.1.g #A single q parameter for everyone
      status.matrix.1$TimeSymptomOnset[first]<-current.time+incubation.period(pathogen=pathogen.1)
      if(runif(1)<bc.1){
      homequarantine.day.1[first]<-status.matrix.1$TimeSymptomOnset[first]
      }
      status.matrix.1$severity[first]<-1
      time.events<-rbind(time.events,c(current.time,1.1,first))
    #}
    #else{
    #  transmission.parameters$q1h[first]<-inf.path.1.h*alpha.as.1 #A single q parameter for everyone
    #  transmission.parameters$q1g[first]<-inf.path.1.g*alpha.as.1 #A single q parameter for everyone
    #  status.matrix.1$severity[first]<-2
    #  time.events<-rbind(time.events,c(current.time,1.2,first))
    #}
    infectives[first]<-1
    contact.time.within$pr.ctc[first]<-ifelse(transmission.parameters$contact_rate_within[first]!=0,rexp(1,transmission.parameters$contact_rate_within[first])+current.time,Inf)       # I generate the next interarrival time for individual i
    contact.time.between$pr.ctc[first]<-rexp(1,transmission.parameters$contact_rate_between[first])+current.time # I generate the next interarrival time for individual i
  }
  
  proposed.individual<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  recovered<-0
  err<-0
  Rt1<-matrix(data = NA, nrow = 1, ncol = 2)
  Rt2<-matrix(data = NA, nrow = 1, ncol = 2)
  
  #When only the first pathogen is present
  while((sum(infectives)>0 & current.time<t.stop) | current.time<t2){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contac
    for (i in which(index.contact.within==1) ){ # for all the individuals that has to propose a global contact
      contact.time.within$pr.ctc[i]<-ifelse(transmission.parameters$contact_rate_within[i]!=0,rexp(1,transmission.parameters$contact_rate_within[i])+current.time,Inf)# I generate the next interarrival time for individual i
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
      
      
      
      #Infection with pathogen 1
#      if (status.matrix.1$infected[infector]==1 & (status.matrix.1$infected[infectee]==0 | (status.matrix.1$infected[infectee]==-1 & reinf==1))){
      if (status.matrix.1$infected[infector]==1 & status.matrix.1$infected[infectee]==0){
        # compute short interaction terms for pathogen.1 (having pathogen 2)
        if (status.matrix.2$infected[infectee]==1){
          short.inter<-sigma12
        }else{
          short.inter<-1
        }
        
        
        #long.inter<-LLImmlev(pathogen.v1 = pathogen.1, pathogen.v2 = pathogen.2, status.matrix.v1 = status.matrix.1, status.matrix.v2 =status.matrix.2, infectee = infectee,lli = lli.2, current.time = current.time, reinf = reinf, typeIC = typeIC, het.vac = het.vac)
        long.inter<-LLImmlev.basic(status.matrix.v2 = status.matrix.2,infectee = infectee, lli = lli.2, current.time = current.time,typeIC = typeIC, t.imm.lim = t.imm.lim, pathogen1 = pathogen.1, pathogen2=pathogen.2 )
        
        # compute long interaction terms for pathogen.1 (recovered from pathogen 2 and/or already contracted pathogen 1)
       # if (status.matrix.2$infected[infectee]==-1 & short.inter==1){
          # long.inter<-long.inter.term.2(t=current.time,status.matrix = status.matrix.2,infectee = infectee)
      #    long.inter<-LLImmlev(pathogen.v1 = pathogen.1, pathogen.v2 = pathogen.2, status.matrix.v1 = status.matrix.1, status.matrix.v2 =status.matrix.2, infectee = infectee,lli = lli.2, current.time = current.time, reinf = reinf, typeIC = typeIC)
      #  }else{
      #    long.inter<-1
      #  }
        
        
        #re-infection term

        #re.inf<-Immlev.1(t=current.time,status.matrix = status.matrix.1,infectee = infectee,pathogen = pathogen.1)
        ifelse(ctc=="g",q<-transmission.parameters$q1g[infector],q<-transmission.parameters$q1h[infector])
        acc.rate.1<-InfMeasure(t= current.time- status.matrix.1$time.of.infection[infector] ,pathogen = pathogen.1)*short.inter*long.inter*q
        if ((ctc=="g" & homequarantine[infectee]==1) | status.matrix.1$infected[infectee]==1){acc.rate.1<-0}
        if (acc.rate.1>1){err<-err+1}
        if (runif(1)<acc.rate.1){
          status.matrix.1$infected[infectee] <- 1 
          status.matrix.1$time.of.infection[infectee] <- current.time
          status.matrix.1$infector[infectee] <- infector
          status.matrix.1$Recovery[infectee]<-current.time+infectious.period.length(pathogen=pathogen.1)
          if (runif(1)<rho.1){ #if symptomatic
            transmission.parameters$q1h[infectee]<-inf.path.1.h #A single q parameter for everyone
            transmission.parameters$q1g[infectee]<-inf.path.1.g #A single q parameter for everyone
            status.matrix.1$TimeSymptomOnset[infectee]<-current.time+incubation.period(pathogen=pathogen.1)
            if (runif(1)<bc.1){
              homequarantine.day.1[infectee]<-status.matrix.1$TimeSymptomOnset[infectee]
            }
            status.matrix.1$severity[infectee]<-1
            time.events<-rbind(time.events,c(current.time,1.1,infectee))
          }else{
            transmission.parameters$q1g[infectee]<-inf.path.1.g*alpha.as.1 #A single q parameter for everyone
            transmission.parameters$q1h[infectee]<-inf.path.1.h*alpha.as.1 #A single q parameter for everyone
            status.matrix.1$severity[infectee]<-2
            time.events<-rbind(time.events,c(current.time,1.2,infectee))
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

     #Infection with pathogen 2
      #if (status.matrix.2$infected[infector]==1 & (status.matrix.2$infected[infectee]==0 | (status.matrix.2$infected[infectee]==-1 & reinf==1))){
      if (status.matrix.2$infected[infector]==1 & status.matrix.2$infected[infectee]==0){
        if (status.matrix.1$infected[infectee]==1){
          short.inter<-sigma21
        }else{
          short.inter<-1
        }
       # if (status.matrix.1$infected[infectee]==-1 & short.inter==1){
          #long.inter<-LLImmlev(pathogen.v1 = pathogen.2, pathogen.v2 = pathogen.1, status.matrix.v1 = status.matrix.2, status.matrix.v2 =status.matrix.1, infectee = infectee,lli=lli.1,current.time = current.time, reinf = reinf, typeIC = typeIC, het.vac = het.vac)
          long.inter<-LLImmlev.basic(status.matrix.v2 = status.matrix.1,infectee = infectee, lli = lli.1, current.time = current.time,typeIC = typeIC, t.imm.lim = t.imm.lim, pathogen1 = pathogen.2, pathogen2=pathogen.1)
          
      #  }else{
      #    long.inter<-1
      #  }
        
        #re-infection term
        #re.inf<-Immlev.2(t=current.time,status.matrix = status.matrix.2,infectee = infectee,pathogen = pathogen.2)
        
        ifelse(ctc=="g",q<-transmission.parameters$q2g[infector],q<-transmission.parameters$q2h[infector])
        acc.rate.2<-InfMeasure(t=(current.time-status.matrix.2$time.of.infection[infector]),pathogen = pathogen.2)*short.inter*long.inter*q
        if ((ctc=="g" & homequarantine[infectee]==1)){acc.rate.2<-0}
        if (acc.rate.2>1){err<-err+1}
        if (runif(1)<acc.rate.2){
          status.matrix.2$infected[infectee] <- 1 
          status.matrix.2$time.of.infection[infectee] <- current.time
          status.matrix.2$infector[infectee] <- infector
          status.matrix.2$Recovery[infectee]<-current.time+infectious.period.length(pathogen=pathogen.2)
          if (runif(1)<rho.2){ #if symptomatic
            transmission.parameters$q2h[infectee]<-inf.path.2.h #A single q parameter for everyone
            transmission.parameters$q2g[infectee]<-inf.path.2.g #A single q parameter for everyone
            status.matrix.2$TimeSymptomOnset[infectee]<-current.time+incubation.period(pathogen=pathogen.2)
            if (runif(1)<bc.2){
              homequarantine.day.2[infectee]<-status.matrix.2$TimeSymptomOnset[infectee]
            }
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
      # compute the long and short interaction terms for pathogen.2
    }
    
    if (next.evts=="HomeQuarantine"){
      current.time<-events$HomeQuarantine
      quarantined.individuals<-which(homequarantine.day.overall==current.time)
      for (k in quarantined.individuals){
        if (k != n & k!= 2*n) {
          temp.ind<- k %% n          
        }else{
          if (k==n){
            homequarantine.day.1[n]<-Inf
            stop.quarantine[n]<-status.matrix.1$Recovery[n]
          }else{
            homequarantine.day.2[n]<-Inf
            stop.quarantine[n]<-status.matrix.2$Recovery[n]
          }
          temp.ind<-n
        }
        if (k>n){
          homequarantine.day.2[temp.ind]<-Inf
          stop.quarantine[temp.ind]<-status.matrix.2$Recovery[temp.ind]
        }
        if (k<n){
          homequarantine.day.1[temp.ind]<-Inf
          stop.quarantine[temp.ind]<-status.matrix.1$Recovery[temp.ind]
        }
        if (homequarantine[temp.ind]==1){ #individual is already in quarantine for the other disease
          stop.quarantine[temp.ind]<-max(status.matrix.2$Recovery[temp.ind],status.matrix.1$Recovery[temp.ind])
        }
        homequarantine[temp.ind]<-1
        contact.time.between$pr.ctc[temp.ind]<-NA
        contact.time.within$pr.ctc[temp.ind]<-NA
        index.contact.between[temp.ind]<-0
        index.contact.within[temp.ind]<-1
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
            Rt2<-comp.RT(status.matrix = status.matrix.2,individual = recovered,Rt=Rt2)
            status.matrix.2$infected[recovered]<--1
            status.matrix.2$Recovery[recovered]<-Inf
            transmission.parameters$contact_rate_within[recovered]<-length(get.neighborhood(HH.network,recovered))
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
            Rt1<-comp.RT(status.matrix = status.matrix.1,individual = recovered,Rt=Rt1)
            status.matrix.1$infected[recovered]<--1
            status.matrix.1$Recovery[recovered]<-Inf
            transmission.parameters$contact_rate_within[recovered]<-length(get.neighborhood(HH.network,recovered))
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
            Rt2<-comp.RT(status.matrix = status.matrix.2,individual = recovered,Rt=Rt2)
            status.matrix.2$infected[recovered]<--1
            status.matrix.2$Recovery[recovered]<-Inf
            transmission.parameters$contact_rate_within[recovered]<-length(get.neighborhood(HH.network,recovered))
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
            Rt1<-comp.RT(status.matrix = status.matrix.1,individual = recovered,Rt=Rt1)
            status.matrix.1$infected[recovered]<--1
            status.matrix.1$Recovery[recovered]<-Inf
            transmission.parameters$contact_rate_within[recovered]<-length(get.neighborhood(HH.network,recovered))
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
          if (runif(1)<bc.2){
            homequarantine.day.2[first]<-status.matrix.2$TimeSymptomOnset[first]
          }
          time.events<-rbind(time.events,c(current.time,2.1,first))
        }else{
          transmission.parameters$q2h[first]<-inf.path.2.h*alpha.as.2 #A single q parameter for everyone
          transmission.parameters$q2g[first]<-inf.path.2.g*alpha.as.2 #A single q parameter for everyone
          status.matrix.2$severity[first]<-2
          time.events<-rbind(time.events,c(current.time,2.2,first))
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
    
    if (next.evts=="NewSeeding1"){
      current.time<-events$NewSeeding1
      events$NewSeeding1<-current.time+t.seed
      not.infected<-which(status.matrix.1$infected!=1)
      if (nSeeds.1<length(not.infected)){
      first.cases<-sample(not.infected,nSeeds.1)
      for (j in first.cases){
        first<-j
        status.matrix.1$infected[first] <- 1 
        status.matrix.1$time.of.infection[first] <- current.time
        status.matrix.1$Recovery[first]<-current.time+infectious.period.length(pathogen=pathogen.1)
        if (runif(1)<rho.1){ #if symptomatic
          transmission.parameters$q1h[first]<-inf.path.1.h #A single q parameter for everyone
          transmission.parameters$q1g[first]<-inf.path.1.g #A single q parameter for everyone
          status.matrix.1$severity[first]<-1
          status.matrix.1$TimeSymptomOnset[first]<-current.time+incubation.period(pathogen=pathogen.2)
          if (runif(1)<bc.1){
            homequarantine.day.1[first]<-status.matrix.1$TimeSymptomOnset[first]
          }
          time.events<-rbind(time.events,c(current.time,1.1,first))
        }else{
          transmission.parameters$q1h[first]<-inf.path.1.h*alpha.as.1 #A single q parameter for everyone
          transmission.parameters$q1g[first]<-inf.path.1.g*alpha.as.1 #A single q parameter for everyone
          status.matrix.1$severity[first]<-2
          time.events<-rbind(time.events,c(current.time,1.2,first))
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
    
    if (next.evts=="NewSeeding2"){
      current.time<-events$NewSeeding1
      events$NewSeeding2<-current.time+t.seed
      not.infected<-which(status.matrix.2$infected!=1)
      if (nSeeds.2<length(not.infected)){
      first.cases<-sample(not.infected,nSeeds.2)
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
          if (runif(1)<bc.2){
            homequarantine.day.2[first]<-status.matrix.2$TimeSymptomOnset[first]
          }
          time.events<-rbind(time.events,c(current.time,2.1,first))
        }else{
          transmission.parameters$q2h[first]<-inf.path.2.h*alpha.as.2 #A single q parameter for everyone
          transmission.parameters$q2g[first]<-inf.path.2.g*alpha.as.2 #A single q parameter for everyone
          status.matrix.2$severity[first]<-2
          time.events<-rbind(time.events,c(current.time,2.2,first))
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
  }
  
  #When also the other pathogen is present.
  time.events<-time.events[-1,]
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  C1<-nSeeds.1
  Y1<-nSeeds.1
  if (t2>0){
  C2<-0
  Y2<-0
  }else{
  C2<-nSeeds.2
  Y2<-nSeeds.2
  }
  last.day<-round(max(time.events[,1]))
  
  for (i in 1:last.day){
    temp.time<-setdiff(which(time.events[,1]>i),which(time.events[,1]>i+1))
    temp.inf.1<-c(which(time.events[temp.time,2]==1.1),which(time.events[temp.time,2]==1.2))
    temp.inf.2<-c(which(time.events[temp.time,2]==2.1),which(time.events[temp.time,2]==2.2))
    temp.time.1<-setdiff(1:length(time.events[,1]),which(time.events[,1]>i+1))
    C1<- c(C1,length((which(time.events[temp.time.1,2]==1.1)))+length((which(time.events[temp.time.1,2]==1.2)))-length((which(time.events[temp.time.1,2]==-1))))
    C2<- c(C2,length((which(time.events[temp.time.1,2]==2.1)))+length((which(time.events[temp.time.1,2]==2.2)))-length((which(time.events[temp.time.1,2]==-2))))
    Y1<- c(Y1,length(temp.inf.1))
    Y2<-c(Y2,length(temp.inf.2))
  }
  
  Fs1<-length(which(time.events[,2]==1.1))+length(which(time.events[,2]==1.2))
  Fs2<-length(which(time.events[,2]==2.1))+length(which(time.events[,2]==2.2))
  
  epi.details<-data.frame("Days"=0:last.day, "Incidence1"=Y1,"Incidence2"=Y2, "Prevalence1"=C1,"Prevalence2"=C2)
  FinalSize<-data.frame("FinalSize1"=Fs1,"FinalSize2"=Fs2)
  PeakIncidence<-data.frame("PeakIncidence1"=max(epi.details$Incidence1),"TimePeakIncidence1"=which(epi.details$Incidence1==max(epi.details$Incidence1))[1],"PeakIncidence2"=max(epi.details$Incidence2),"TimePeakIncidence2"=which(epi.details$Incidence2==max(epi.details$Incidence2))[1] )
  PeakPrevalence<-data.frame("PeakPrevalence1"=max(epi.details$Prevalence1),"TimePeakPrevalence1"=which(epi.details$Prevalence1==max(epi.details$Prevalence1))[1],"PeakPrevalence2"=max(epi.details$Prevalence2),"TimePeakPrevalence2"=which(epi.details$Prevalence2==max(epi.details$Prevalence2))[1] )
  return(list(time.events=time.events, status.matrix.1=status.matrix.1, status.matrix.2=status.matrix.2,epi.details=epi.details, FinalSize=FinalSize, PeakIncidence=PeakIncidence, PeakPrevalence=PeakPrevalence,Rt1=Rt1,Rt2=Rt2))
}














