
############################################################################
############################################################################
### This script contains a set of functions necessary for the simulation ###
### Called by the script 'main.coinfection.new.R'                        ###
############################################################################
############################################################################



#####################################################################
### This function defines the number of days before symptom onset ###
#####################################################################

incubation.period <- function(pathogen){
  if(pathogen == "COVID-19" | pathogen == "DELTA" | pathogen== "OMICRON"){
    #return(rlnorm(1,meanlog = log(5.2), sdlog = log(1.7)))
    return(5.288462)
  }
  if(pathogen == "FLU-A"){
    #return(rlnorm(1,meanlog = log(1.4), sdlog = log(1.51))) 
    return(2)
  }
  if(pathogen == "RSV"){
    return(rlnorm(1,meanlog = log(4.4), sdlog = log(1.24))) 
  }
}


#########################################################################
### This function defines the length of the infectious period in days ###
#########################################################################

infectious.period.length <- function(pathogen){
  if(pathogen == "COVID-19" | pathogen == "OMICRON" | pathogen == "DELTA"){
    return(12)
  }
  if(pathogen == "FLU-A"){
    return(8)  
  }
}


#################################################################################
### This function defines the infectiousness measure.                         ###
### It describes how the probability of transmission varies over the course   ###
### of infection.                                                             ###
### It is estimated by using a viral load curve                               ###
### (see the example for influenza using data reported in Carrat et al. 2008) ###
#################################################################################

InfMeasure <- function(t, pathogen){
  if(pathogen == "COVID-19" | pathogen == "DELTA" | pathogen == "OMICRON"){
    return(dgamma(t,shape = 12, rate = 2.08)/ (pgamma(12,shape = 12,rate = 2.08)))
  }
  if(pathogen == "FLU-A"){
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
  if(pathogen == "RSV"){
    return(dgamma(t,shape = 15, rate = 2.6)/ (pgamma(12,shape = 15,rate = 2.6)) )
  }
}


#################################################################
### This function defines the vaccine effectiveness over time ###
#################################################################

VaccineEffectiveness <- function(t,typeIC){
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



##################################################################
### This function computes the susceptibility of an individual ###
### depending on infectious history                            ###
##################################################################

LLImmlev.basic <- function(status.matrix.v2, infectee, long.int, current.time, typeIC,
                           t.imm.lim, pathogen1, pathogen2){ #pathogen.v1 is the infection the infectee might catch
  value <- 1
  if(status.matrix.v2$infected[infectee] == -1){
    #t.sinc.inf<-current.time-(status.matrix.v2$time.recovery[infectee]+infectious.period.length(pathogen = pathogen2))
    t.sinc.inf <- current.time - (status.matrix.v2$time.infection[infectee] 
                                  + infectious.period.length(pathogen = pathogen2)) # here we are assuming a constant IPL
    if (t.sinc.inf<t.imm.lim){
      if (typeIC==1){
        value<-long.int
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
          value<-long.int
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


#Computing Rt - started but haven't look at this anymore
# Rt
comp.RT<-function(status.matrix,individual,Rt){
  infectees<-which(status.matrix$infected.by==individual)
  Rt.temp<-0
  if (length(infectees)>0){
    for (i1 in 1:length(infectees)){
      if (status.matrix$time.infection[infectees[i1]]>status.matrix$time.infection[individual]){
        Rt.temp<-Rt.temp+1
      }
    }
  }
  Rt<-rbind(Rt,c(status.matrix$time.infection[individual],Rt.temp))
  return(Rt)
}


###############################################################################################
### This function runs the simulation of the IBM.                                           
###                                                                                          
### Input parameters:                                                                       
### t2:                   time at which pathogen 2 is introduced in the population          
### sigma[1]:              short-term interaction parameter: acquiring 2 while having 1      
###                       (if >1 cooperative effect - if <1 competing)                      
### sigma21:              short-term interaction parameter: acquiring 1 while having 2      
###                       (if >1 cooperative effect - if <1 competing) 
### prop.immune:          proportion of immune cases (not used at the moment)
### n.seeds.1:             number of initial cases for path 1
### n.seeds.2:             number of initial cases for path 2
### rho[1]:                probability of being symptomatic for path 1
### rho[2]:                probability of being symptomatic for path 2
### alpha.as[1]:           relative infectiousness of asymptomatic cases (pathogen1)
### alpha.as[2]:           relative infectiousness of asymptomatic cases (pathogen2)
### netw:                 type of household network considered - Synthetic or ERGM
### n.vertex:             number of vertexes 
### n.networks:           number of simulated networks
### R.1:                  reproduction number path 1 (household R*)
### R.2:                  reproduction number path 2 (household R*)
### ratio.qhqg:           ratio transmission probability given household contact 
###                       over global contacts
### long.int[1]:                long-term interaction parameter: acquiring 2 while having  
###                       experienced (and recovered from) 1 
### long.int[2]:                long-term interaction parameter: acquiring 1 while having  
###                       experienced (and recovered from) 2
### pathogen[1]:           character variable identifying pathogen 1
### pathogen[2]:           character variable identifying pathogen 2
### contact.reduction:    parameter multiplying the household contact rate after home 
###                       isolation
### t.stop:               time at which simulations stop
### t.seed:               time of additional seeding
### behavior.change.1:                 proportion of individuals changing behavior (home isolation) 
###                       after being infected with pathogen 1
### behavior.change.2:                 proportion of individuals changing behavior (home isolation) 
###                       after being infected with pathogen 2
### reinfection:                boolean identifying whether someone can be re-infected with 
###                       the same pathogen (1 yes, 0 no)
### typeIC:               ID for different type of waning of immunity
### contact.reduction.TP: contact reduction value set to identify transmission rates 
###                       (household and global) linked to a specific R*
### behavior.change.1.TP:              behavior change value (for pathogen 1) set to identify transmission 
###                       rates (household and global) linked to a specific R*
### behavior.change.2.TP:              behavior change value (for pathogen 2) set to identify transmission 
###                       rates (household and global) linked to a specific R*
### het.vac:              boolean for heterologous effects (1 yes 0 no) - Not used currently
### t.imm.lim:            parameter to define the length of immunity that have the same overall 
###                       "effect" (area underneath the curve)
### decrease.gc:               decrease in the  number of global contact rates compared to baseline

### Output parameters:
### time events: three-column matrix identifying one of the following events(second column):
###                 1.1 and infection with pathogen 1 occurs
###                 1.2 and infection with pathogen 2 occurs
###                 0.1 an individual recovers from an infection with pathogen 1
###                 0.2 an individual recovers from an infection with pathogen 2
###                     the time at which the event occur (first column) 
###                     and the individual (third column)
### status matrix.1/status matrix.1: matrix that represents in each row and individual, and 
###                                  infection related characteristics 
###############################################################################################

sim.multipathogen <- function(HH.network, t2, lambda.g, sigma, prop.immune, n.seeds, 
                              rho, inf.path.h, inf.path.g, 
                              alpha.as, long.int, pathogen, contact.reduction,
                              t.stop, t.seed, behavior.change, reinfection, typeIC, het.vac, t.imm.lim){
  
  n <- network.size(HH.network)
  hh.id <- HH.network %v% "hh_id"
  hh.size<- HH.network %v% "hh_size"
  
  
  ### Create data structures to keep track of the results
  status.matrix <- data.frame(infected = rep(FALSE,n),
                              time.infection = NA,
                              infected.by = NA,
                              severity = 0, # 1 symptomatic, 2 asymptomatic
                              time.symptom.onset = Inf,
                              immunity = 0, # 0 no immunity, 1 vaccinated
                              time.recovery = Inf)
  status = list(status.matrix, status.matrix)
  
  time.recovery.vector.1 <- time.recovery.vector.2 <- rep(Inf,n) #vector giving the time.recovery times

  # dataframe containing the next time of each possible event types
  events <- data.frame(contact = Inf, 
                       home.quarantine = Inf, 
                       recovery = Inf, 
                       new.pathogen = t2, 
                       new.seeding.1 = t.seed,
                       new.seeding.2 = t2 + t.seed)
  # vector indicating who is infectious at the current time: 1 infectious 0 non infectious
  infectives <- rep(FALSE,n) 
  # vector that selects the individuals that have to propose a new social contact(local) - 1 yes 0 no
  index.contact.h <- rep(FALSE,n) 
  # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  index.contact.g <- rep(FALSE,n) 
  # ?
  time.events = data.frame(time = NA,
                           type = NA,
                           id = NA)  
  # Set the current time
  current.time <- 0
  
  # transmission parameter dataframe: each line is an individual, the first column is the ID, the second and third the transmisssion parameter given  household or global contacts for pathogen 1, and fourth and fifth for pathogen 2
  # 6th, 7th columns the contact rates for within and between household contacts, and last a vector checking whether a person is susceptible
  transmission.parameters <- data.frame(id = 1:n, 
                                        q.1.h = NA, 
                                        q.1.g = NA, 
                                        q.2.h = NA, 
                                        q.2.g = NA, 
                                        contact.rate.h = sapply(1:n, function(j) length(get.neighborhood(HH.network, j))),
                                        contact.rate.g = lambda.g, 
                                        susceptibility = 1)  
  
  #Proportion of immune - To be adapted
  
  # if (prop.immune>0){
  #   if (pathogen[1]=="DELTA" & pathogen[2]=="OMICRON"){
  #     immuned.individuals<-sample(1:n,round(prop.immune*n))
  #     status.matrix.1$immunity[immuned.individuals]<-1
  #     status.matrix.2$immunity<-status.matrix.1$immunity
  #   }else{
  #     immuned.individuals<-sample(1:n,round(prop.immune*n))
  #     status.matrix.1$immunity[immuned.individuals]<-1
  #     immuned.individuals<-sample(1:n,round(prop.immune*n))
  #     status.matrix.2$immunity[immuned.individuals]<-1
  #   }
  # }
  
  # individuals can home isolate when they show symptoms - if so, they stay in iso till the end of the infectious period (for sake of simplicity)
  
  home.quarantine = data.frame(id = 1:n,
                          quarantine = FALSE,
                          start = Inf,
                          end = Inf)
  home.quarantine = list(home.quarantine, home.quarantine)
  
  # matrix containing the proposed time of the next contact (second column) and the contact individual (third column) for within household contacts
  
  next.contact = list(h = data.frame(id = 1:n, 
                                     time.next.contact = Inf, 
                                     person.time.next.contact = NA),
                      g = data.frame(id = 1:n, 
                                     time.next.contact = Inf, 
                                     person.time.next.contact = NA))
  
  ##################################################################
  ### STEP 1
  ### Introduce the first infections
  ### The individuals are randomly chosen in the population (among susceptibles).
  ### We assume that epidemics always start in household of size 2 (this is a random choice that decreases stochasticity).
  ### We assume that index cases are always symptomatic.
  first.cases.1 <- identify.initial.cases(pathogen[1])
  
  for(j in first.cases.1){
    initial.infection(id = j, path = pathogen[1])
  }
  
  ### ? 
  proposed.individual <- 0
  temp.contact.time <- 0
  indiv.prop.ctc <- 0
  recovered <- 0
  err <- 0
  Rt1 <- matrix(data = NA, nrow = 1, ncol = 2)
  Rt2 <- matrix(data = NA, nrow = 1, ncol = 2)
  ### ?
  
  ##################################################################
  ### STEP 2
  ### Run the epidemic  
  ### Continue the epidemic as long as there are infected individuals or pathogen 2 still has to be introduced and as long
  ### as t.stop has not been reached.

  #while((sum(infectives)>0 & current.time<t.stop) | current.time<t2){ #while there are still infectives, we are within the t.stop
  while((sum(status[[1]]$infected + status[[2]]$infected) > 0 & current.time < t.stop) | current.time < t2){
    
    ### Phase 1: individuals that has to, propose a new social contact
    generate.next.contact()
    
    ###  Phase 2: identify the next event: select the minimum time among the events that can occur
    next.event = identify.next.event()
    
    ### PHASE 3: the event happens
    
    
    if (next.event=="contact"){
      current.time<-events$contact
      if (length(min(contact.time.overall, na.rm = T))>1){ #when two contacts happen at the same time select one of the two
        selected.ctc<-sample(which(contact.time.overall==current.time),1) 
        if (selected.ctc!=n & selected.ctc!=2*n){
          infected.by<- selected.ctc %% n
        }else{
          infected.by<- n
        }
        if (selected.ctc<=n){
          infectee.pool<-get.neighborhood(HH.network,infected.by)
          if (length(infectee.pool)>1){
            infectee<-sample(infectee.pool,1) #pick a random individual in the class
          }
          if(length(infectee.pool)==1){
            infectee<-infectee.pool #pick a random individual in the class
          }
          if(length(infectee.pool)==0){
            infectee<-infected.by #just a trick to have acceptance rate 0 (infected.by is not susceptible)
          }
          index.contact.h[infected.by]<-1
          next.contact.h$time.next.contact[infected.by]<-NA
          ctc<-"hh"
        }else{
          infected.by<-which(next.contact.g$time.next.contact == events$contact) 
          hh.members.temp<-which(hh.id==hh.id[infected.by])
          infectee.pool<- setdiff(1:n,hh.members.temp) #individuals not in the same class
          infectee<-sample(infectee.pool,1) #pick a random individual not in the class (and not a teacher)
          index.contact.g[infected.by]<-1
          next.contact.g$time.next.contact[infected.by]<-NA
          ctc<-"g"
        }
      }else{
        if (length(which(events$contact==next.contact.h$time.next.contact))>0){ #if it is a within contact
          infected.by<-which(next.contact.h$time.next.contact ==events$contact)
          infectee.pool<-get.neighborhood(HH.network,infected.by) # the pool of susceptible is the household
          if (length(infectee.pool)>1){
            infectee<-sample(infectee.pool,1) #pick a random individual in the class
          }
          if(length(infectee.pool)==1){
            infectee<-infectee.pool #pick a random individual in the class
          }
          if(length(infectee.pool)==0){
            infectee<-infected.by #just a trick to have acceptance rate 0 (infected.by is not susceptible)
          }
          index.contact.h[infected.by]<-1
          next.contact.h$time.next.contact[infected.by]<-NA
          ctc<-"hh"
        }else{
          infected.by<-which(next.contact.g$time.next.contact == events$contact) 
          hh.members.temp<-which(hh.id==hh.id[infected.by])
          infectee.pool<- setdiff(1:n,hh.members.temp) #individuals not in the same class
          infectee<-sample(infectee.pool,1) #pick a random individual not in the class (and not a teacher)
          index.contact.g[infected.by]<-1
          next.contact.g$time.next.contact[infected.by]<-NA
          ctc<-"g"
        }
      }
      
      
      
      #Infection with pathogen 1
      if (status.matrix.1$infected[infected.by]== TRUE & status.matrix.1$infected[infectee]== FALSE){
        # compute short interaction terms for pathogen.1 (having pathogen 2)
        if (status.matrix.2$infected[infectee]==TRUE){
          short.inter<-sigma[1]
        }else{
          short.inter<-1
        }
        # compute long-term interactions (i.e., the potential infectee already experienced the other infection, or that inefction)
        long.inter<-LLImmlev.basic(status.matrix.v2 = status.matrix.2,infectee = infectee, long.int = long.int[2], current.time = current.time,typeIC = typeIC, t.imm.lim = t.imm.lim, pathogen1 = pathogen[1], pathogen2=pathogen[2] )
        # select the transmission probability related to global or local contacts
        ifelse(ctc=="g",q<-transmission.parameters$q.1.g[infected.by],q<-transmission.parameters$q.1.h[infected.by])
        #Acceptance rate is the probability that there will be the infection. This is composed by q,short and long term interaction, and the infectiousness measure that describe how likely is that it will happen in that moment
        acc.rate.1<-InfMeasure(t= current.time- status.matrix.1$time.infection[infected.by] ,pathogen = pathogen[1])*short.inter*long.inter*q
        if ((ctc=="g" & home.quarantine[infectee]==1) | status.matrix.1$infected[infectee]==TRUE){acc.rate.1<-0} # if the contacted person is in home quarantine there is no contact 
        if (acc.rate.1>1){err<-err+1} #keep track that the acc.rate does not exceed 1 - since it is a prob value this would not make much sense
        if (runif(1)<acc.rate.1){ # if infection, adapt the status matrix
          status.matrix.1$infected[infectee] <- TRUE
          status.matrix.1$time.infection[infectee] <- current.time
          status.matrix.1$infected.by[infectee] <- infected.by
          status.matrix.1$time.recovery[infectee]<-current.time+infectious.period.length(pathogen=pathogen[1])
          if (runif(1)<rho[1]){ #if symptomatic
            transmission.parameters$q.1.h[infectee]<-inf.h[1] #A single q parameter for everyone
            transmission.parameters$q.1.g[infectee]<-inf.g[1] #A single q parameter for everyone
            status.matrix.1$time.symptom.onset[infectee]<-current.time+incubation.period(pathogen=pathogen[1])
            if (runif(1)<behavior.change.1){
              home.quarantine.start.1[infectee]<-status.matrix.1$time.symptom.onset[infectee]
            }
            status.matrix.1$severity[infectee]<-1
            time.events<-rbind(time.events,c(current.time,1.1,infectee))
          }else{
            transmission.parameters$q.1.g[infectee]<-inf.g[1]*alpha.as[1] #A single q parameter for everyone
            transmission.parameters$q.1.h[infectee]<-inf.h[1]*alpha.as[1] #A single q parameter for everyone
            status.matrix.1$severity[infectee]<-2
            time.events<-rbind(time.events,c(current.time,1.2,infectee))
          }
          if (infectives[infectee]==0){
            infectives[infectee]<-1
            next.contact.h$time.next.contact[infectee]<-ifelse(transmission.parameters$contact.rate.h[infectee]!=0,rexp(1,transmission.parameters$contact.rate.h[infectee])+current.time,Inf)       # I generate the next interarrival time for individual i
            if (home.quarantine[infectee]==0){
              next.contact.g$time.next.contact[infectee]<-rexp(1,transmission.parameters$contact.rate.g[infectee])+current.time # I generate the next interarrival time for individual i
            }
          }
        }
        
      }

     #Infection with pathogen 2 - similar to above. To note: if the infected.by is infected with both, both these probability will be computed and potentially both infections can be transmitted
      if (status.matrix.2$infected[infected.by]== TRUE & status.matrix.2$infected[infectee]==FALSE){
        if (status.matrix.1$infected[infectee]==TRUE){
          short.inter<-sigma21
        }else{
          short.inter<-1
        }
          long.inter<-LLImmlev.basic(status.matrix.v2 = status.matrix.1,infectee = infectee, long.int = long.int[1], current.time = current.time,typeIC = typeIC, t.imm.lim = t.imm.lim, pathogen1 = pathogen[2], pathogen2=pathogen[1])
          
        ifelse(ctc=="g",q<-transmission.parameters$q.2.g[infected.by],q<-transmission.parameters$q.2.h[infected.by])
        acc.rate.2<-InfMeasure(t=(current.time-status.matrix.2$time.infection[infected.by]),pathogen = pathogen[2])*short.inter*long.inter*q
        if ((ctc=="g" & home.quarantine[infectee]==1)){acc.rate.2<-0}
        if (acc.rate.2>1){err<-err+1}
        if (runif(1)<acc.rate.2){
          status.matrix.2$infected[infectee] <- TRUE 
          status.matrix.2$time.infection[infectee] <- current.time
          status.matrix.2$infected.by[infectee] <- infected.by
          status.matrix.2$time.recovery[infectee]<-current.time+infectious.period.length(pathogen=pathogen[2])
          if (runif(1)<rho[2]){ #if symptomatic
            transmission.parameters$q.2.h[infectee]<-inf.h[2] #A single q parameter for everyone
            transmission.parameters$q.2.g[infectee]<-inf.g[2] #A single q parameter for everyone
            status.matrix.2$time.symptom.onset[infectee]<-current.time+incubation.period(pathogen=pathogen[2])
            if (runif(1)<behavior.change.2){
              home.quarantine.start.2[infectee]<-status.matrix.2$time.symptom.onset[infectee]
            }
            status.matrix.2$severity[infectee]<-1
            time.events<-rbind(time.events,c(current.time,2.1,infectee))
          }else{
            transmission.parameters$q.2.g[infectee]<-inf.g[2]*alpha.as[2] #A single q parameter for everyone
            transmission.parameters$q.2.h[infectee]<-inf.h[2]*alpha.as[2] #A single q parameter for everyone
            status.matrix.2$severity[infectee]<-2
            time.events<-rbind(time.events,c(current.time,2.2,infectee))
          }
          if (infectives[infectee]==0){
            infectives[infectee]<-1
            next.contact.h$time.next.contact[infectee]<-ifelse(transmission.parameters$contact.rate.h[infectee]!=0,rexp(1,transmission.parameters$contact.rate.h[infectee])+current.time,Inf)       # I generate the next interarrival time for individual i
            if (home.quarantine[infectee]==0){
              next.contact.g$time.next.contact[infectee]<-rexp(1,transmission.parameters$contact.rate.g[infectee])+current.time # I generate the next interarrival time for individual i
            }
          }
        }
        
      }
    }
    
    if (next.event=="home.quarantine"){ #next event is home isolation
      current.time<-events$home.quarantine
      quarantined.individuals<-which(home.quarantine.start.overall==current.time)
      for (k in quarantined.individuals){ # we need to understand whether is an individual sick by pathogen 1 or 2
        if (k != n & k!= 2*n) {
          temp.ind<- k %% n          
        }else{
          if (k==n){
            home.quarantine.start.1[n]<-Inf
            stop.quarantine[n]<-status.matrix.1$time.recovery[n] # home isolated till time.recovery
          }else{
            home.quarantine.start.2[n]<-Inf
            stop.quarantine[n]<-status.matrix.2$time.recovery[n]
          }
          temp.ind<-n
        }
        if (k>n){
          home.quarantine.start.2[temp.ind]<-Inf
          stop.quarantine[temp.ind]<-status.matrix.2$time.recovery[temp.ind]
        }
        if (k<n){
          home.quarantine.start.1[temp.ind]<-Inf
          stop.quarantine[temp.ind]<-status.matrix.1$time.recovery[temp.ind]
        }
        if (home.quarantine[temp.ind]==1){ #individual is already in quarantine for the other disease
          stop.quarantine[temp.ind]<-max(status.matrix.2$time.recovery[temp.ind],status.matrix.1$time.recovery[temp.ind])
        }
        home.quarantine[temp.ind]<-1
        next.contact.g$time.next.contact[temp.ind]<-NA
        next.contact.h$time.next.contact[temp.ind]<-NA
        index.contact.g[temp.ind]<-0
        index.contact.h[temp.ind]<-1
        transmission.parameters$contact.rate.h[temp.ind]<-transmission.parameters$contact.rate.h[temp.ind]*contact.reduction
      }
    }
    
    if (next.event=="time.recovery"){ #time.recovery from infection
      current.time<-events$recovery
      temp.recovered<-which(time.recovery.vector.overall==events$recovery)
      for (recovered in temp.recovered){
        if (recovered!= n & recovered!=n*2){ # whether recovering from pathogen 1 or 2
          if (recovered > n){
            recovered<- recovered %% n
            Rt2<-comp.RT(status.matrix = status.matrix.2,individual = recovered,Rt=Rt2)
            status.matrix.2$infected[recovered]<-- TRUE
            status.matrix.2$time.recovery[recovered]<-Inf
            transmission.parameters$contact.rate.h[recovered]<-length(get.neighborhood(HH.network,recovered))
            time.events<-rbind(time.events,c(current.time,-2,recovered))
            if (status.matrix.1$infected[recovered]!=TRUE){
              infectives[recovered]<-0
              next.contact.g$time.next.contact[recovered]<-NA #be sure to reomove all the contacts that these person will make after time.recovery (if time.recovery comes before nextcontact)
              next.contact.h$time.next.contact[recovered]<-NA 
              index.contact.h[recovered]<-0 #this person is not making any other contact (this will happen only if infected agai)
              index.contact.g[recovered]<-0
            }else{
              if (home.quarantine[recovered]==0){ #remove from home isolation - so he/she can be reached by global contacts
                index.contact.g[recovered]<-1
              }
            }
          }else{
            Rt1<-comp.RT(status.matrix = status.matrix.1,individual = recovered,Rt=Rt1)
            status.matrix.1$infected[recovered]<--TRUE
            status.matrix.1$time.recovery[recovered]<-Inf
            transmission.parameters$contact.rate.h[recovered]<-length(get.neighborhood(HH.network,recovered))
            time.events<-rbind(time.events,c(current.time,-1,recovered))
            if (status.matrix.2$infected[recovered]!=TRUE){
              infectives[recovered]<-0
              next.contact.g$time.next.contact[recovered]<-NA
              next.contact.h$time.next.contact[recovered]<-NA
              index.contact.h[recovered]<-0
              index.contact.g[recovered]<-0
            }else{
              if (home.quarantine[recovered]==0){
                index.contact.g[recovered]<-1
              }
            }
          }
        }else{
          if (recovered == 2*n){
            recovered<- n
            Rt2<-comp.RT(status.matrix = status.matrix.2,individual = recovered,Rt=Rt2)
            status.matrix.2$infected[recovered]<--TRUE
            status.matrix.2$time.recovery[recovered]<-Inf
            transmission.parameters$contact.rate.h[recovered]<-length(get.neighborhood(HH.network,recovered))
            time.events<-rbind(time.events,c(current.time,-2,recovered))
            if (status.matrix.1$infected[recovered]!=TRUE){
              infectives[recovered]<-0
              next.contact.g$time.next.contact[recovered]<-NA
              next.contact.h$time.next.contact[recovered]<-NA
              index.contact.h[recovered]<-0
              index.contact.g[recovered]<-0
            }else{
              if (home.quarantine[recovered]==0){
                index.contact.g[recovered]<-1
              }
            }
          }else{
            Rt1<-comp.RT(status.matrix = status.matrix.1,individual = recovered,Rt=Rt1)
            status.matrix.1$infected[recovered]<--TRUE
            status.matrix.1$time.recovery[recovered]<-Inf
            transmission.parameters$contact.rate.h[recovered]<-length(get.neighborhood(HH.network,recovered))
            time.events<-rbind(time.events,c(current.time,-1,recovered))
            if (status.matrix.2$infected[recovered]!=TRUE){
              infectives[recovered]<-0
              next.contact.g$time.next.contact[recovered]<-NA
              next.contact.h$time.next.contact[recovered]<-NA
              index.contact.h[recovered]<-0
              index.contact.g[recovered]<-0
            }else{
              if (home.quarantine[recovered]==0){
                index.contact.g[recovered]<-1
              }
            }
          }
        }
        if (stop.quarantine[recovered]==current.time){
          home.quarantine[recovered]<-0
          stop.quarantine[recovered]<-Inf
        }
      }
    }
    
    if(next.event == "new.pathogen"){ 
      # update the current time: 
      current.time <- events$new.pathogen
      # identify the individuals to be infected
      first.cases.2 <- identify.initial.cases(pathogen[2])
      # infect thpse individuals
      for(j in first.cases.2){
        initial.infection(id = j, path = pathogen[2])
      }
      # the event has been introduced, the time needs to be set to Inf. 
      events$new.pathogen <- Inf
    }
    
    
    
    
    if(next.event == "new.seeding.1"){
      current.time <- events$new.seeding.1
      events$new.seeding.1<-current.time+t.seed
      not.infected<-which(status.matrix.1$infected!=TRUE)
      if (n.seeds.1<length(not.infected)){
      first.cases<-sample(not.infected,n.seeds.1)
      for (j in first.cases){
        first<-j
        status.matrix.1$infected[first] <- TRUE
        status.matrix.1$time.infection[first] <- current.time
        status.matrix.1$time.recovery[first]<-current.time+infectious.period.length(pathogen=pathogen[1])
        if (runif(1)<rho[1]){ #if symptomatic
          transmission.parameters$q.1.h[first]<-inf.h[1] #A single q parameter for everyone
          transmission.parameters$q.1.g[first]<-inf.g[1] #A single q parameter for everyone
          status.matrix.1$severity[first]<-1
          status.matrix.1$time.symptom.onset[first]<-current.time+incubation.period(pathogen=pathogen[2])
          if (runif(1)<behavior.change.1){
            home.quarantine.start.1[first]<-status.matrix.1$time.symptom.onset[first]
          }
          time.events<-rbind(time.events,c(current.time,1.1,first))
        }else{
          transmission.parameters$q.1.h[first]<-inf.h[1]*alpha.as[1] #A single q parameter for everyone
          transmission.parameters$q.1.g[first]<-inf.g[1]*alpha.as[1] #A single q parameter for everyone
          status.matrix.1$severity[first]<-2
          time.events<-rbind(time.events,c(current.time,1.2,first))
        }
        if (infectives[first]==0){
          infectives[first]<-1
          next.contact.h$time.next.contact[first]<-ifelse(transmission.parameters$contact.rate.h[first]!=0,rexp(1,transmission.parameters$contact.rate.h[first])+current.time,Inf)       # I generate the next interarrival time for individual i
          if (home.quarantine[first]==0){
            next.contact.g$time.next.contact[first]<-rexp(1,transmission.parameters$contact.rate.g[first])+current.time # I generate the next interarrival time for individual i
          }
        }
      }
      }
    }
    
    if (next.event=="new.seeding.2"){
      current.time<-events$new.seeding.1
      events$new.seeding.2<-current.time+t.seed
      not.infected<-which(status.matrix.2$infected!=TRUE)
      if (n.seeds.2<length(not.infected)){
      first.cases<-sample(not.infected,n.seeds.2)
      for (j in first.cases){
        first<-j
        status.matrix.2$infected[first] <- TRUE 
        status.matrix.2$time.infection[first] <- current.time
        status.matrix.2$time.recovery[first]<-current.time+infectious.period.length(pathogen=pathogen[2])
        if (runif(1)<rho[2]){ #if symptomatic
          transmission.parameters$q.2.h[first]<-inf.h[2] #A single q parameter for everyone
          transmission.parameters$q.2.g[first]<-inf.g[2] #A single q parameter for everyone
          status.matrix.2$severity[first]<-1
          status.matrix.2$time.symptom.onset[first]<-current.time+incubation.period(pathogen=pathogen[2])
          if (runif(1)<behavior.change.2){
            home.quarantine.start.2[first]<-status.matrix.2$time.symptom.onset[first]
          }
          time.events<-rbind(time.events,c(current.time,2.1,first))
        }else{
          transmission.parameters$q.2.h[first]<-inf.h[2]*alpha.as[2] #A single q parameter for everyone
          transmission.parameters$q.2.g[first]<-inf.g[2]*alpha.as[2] #A single q parameter for everyone
          status.matrix.2$severity[first]<-2
          time.events<-rbind(time.events,c(current.time,2.2,first))
        }
        if (infectives[first]==0){
          infectives[first]<-1
          next.contact.h$time.next.contact[first]<-ifelse(transmission.parameters$contact.rate.h[first]!=0,rexp(1,transmission.parameters$contact.rate.h[first])+current.time,Inf)       # I generate the next interarrival time for individual i
          if (home.quarantine[first]==0){
            next.contact.g$time.next.contact[first]<-rexp(1,transmission.parameters$contact.rate.g[first])+current.time # I generate the next interarrival time for individual i
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
  
  # compute some summary measures that will be given as output
  
  C1<-n.seeds.1
  Y1<-n.seeds.1
  if (t2>0){
  C2<-0
  Y2<-0
  }else{
  C2<-n.seeds.2
  Y2<-n.seeds.2
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

##########################################################################
### This function choses the initial infections for a certain pathogen ###
##########################################################################

identify.initial.cases = function(path){
  p <- ifelse(path == pathogen[1], 1, 2)
  potential.seeds <- which(hh.size == 2)
  first.cases <<- sample(potential.seeds[which(status[[p]][potential.seeds, "infected"] == FALSE)], n.seeds[p])
  return(first.cases)
}

#########################################################################
### This function does what has to happen during an event "infection" ###
#########################################################################

initial.infection = function(id, path){
  # index for the pathogen
  p = ifelse(path == pathogen[1], 1, 2)
  # update the status
  status[[p]]$infected[id] <<- TRUE
  status[[p]]$time.infection[id] <<- current.time
  status[[p]]$time.recovery[id] <<- current.time + infectious.period.length(pathogen = pathogen[p])
  if(runif(1) < rho[p]){ # if the individual is symptomatic
    if(p == 1){
      transmission.parameters$q.1.h[id] <<- inf.h[1]
      transmission.parameters$q.1.g[id] <<- inf.g[1]
    }else{
      transmission.parameters$q.2.h[id] <<- inf.h[2]
      transmission.parameters$q.2.g[id] <<- inf.g[2]
    }
    status[[p]]$severity[id] = 1
    status[[p]]$time.symptom.onset[id] <<- current.time+incubation.period(pathogen=pathogen[p])
    if(runif(1) < behavior.change[p]){
      home.quarantine[[p]]$start[id] <<- status[[p]]$time.symptom.onset[id]
    }
    time.events <<- rbind(time.events, c(current.time, "symptomaitc infection", id))
  }else{
    if(p == 1){
      transmission.parameters$q.1.h[id] <<- inf.h[1] * alpha.as[1]
      transmission.parameters$q.1.g[id] <<- inf.g[1] * alpha.as[1]
    }else{
      transmission.parameters$q.2.h[id] <<- inf.h[2] * alpha.as[2]
      transmission.parameters$q.2.g[id] <<- inf.g[2] * alpha.as[2]
    }
    status[[p]]$severity[id] = 2
    time.events <<- rbind(time.events, c(current.time, "asymptomaitc infection", id))
    
  }
  # generate the time of the next household contact
  next.contact$h$time.next.contact[id] <<- ifelse(transmission.parameters$contact.rate.h[id] != 0, # check whether there are any household members
                                                 rexp(1, transmission.parameters$contact.rate.h[id]) + current.time,
                                                 Inf)
  # if the individual is not in home quarantine, generate the time of the next global contact
  if(home.quarantine[[1]]$quarantine[id] == FALSE & home.quarantine[[2]]$quarantine[id] == FALSE){
    next.contact$g$time.next.contact[id] <<- rexp(1, transmission.parameters$contact.rate.g[id]) + current.time
  }
  infectives[id] <<- TRUE
}

#################################################
### This function generates the next contacts ###
#################################################

generate.next.contact = function(){
  # For each individual that has to propose a household contact, generate the time of the next contact
  for(j in which(index.contact.h == TRUE)){ 
    next.contact$h$time.next.contact[j] <<- ifelse(transmission.parameters$contact.rate.h[j] != 0, # check whether there are any household members
                                                   rexp(1, transmission.parameters$contact.rate.h[j]) + current.time,
                                                   Inf)
    # this individual no longer has to propose a new contact
    index.contact.h[j] <<- FALSE
  }
  
  # For each individual that has to propose a global contact, generate the time of the next contact
  for(j in which(index.contact.g == TRUE)){
    next.contact$g$time.next.contact[j] <<- rexp(1, transmission.parameters$contact.rate.g[j]) + current.time
    # this individual no longer has to propose a new contact
    index.contact.g[j] <<- FALSE
  }
}

###############################################
### This function identifies the next event ###
###############################################

identify.next.event = function(){
  # identify the next contact among all household and global contacts
  events$contact <<- min(c(next.contact$h$time.next.contact, next.contact$g$time.next.contact), na.rm = TRUE)
  # identify the next home.quarantine time for both pathogens
  events$home.quarantine <<- min(c(home.quarantine[[1]]$start, home.quarantine[[2]]$start), na.rm = TRUE)
  # identify the next recovery time for both pathogens
  events$recovery <<- min(c(status[[1]]$time.recovery, status[[2]]$time.recovery), na.rm = TRUE)
  
  # select the next event(s). If more than one event: sample one at random
  next.event <- colnames(events)[which(min(events) == events)]
  if (length(next.event) > 1){
    next.event <- sample(next.events ,1)
  }
  return(next.event)
}












