
############################################################################
############################################################################
### This script contains a set of functions necessary for the simulation ###
### Called by the script 'main.coinfection.new.R'                        ###
############################################################################
############################################################################

###############################################################################################
### This function runs the simulation of the IBM.                                           
###                                                                                          
### Input parameters:                                                                       
### t2:                   time at which pathogen 2 is introduced in the population          
### sigma[1]:             short-term interaction parameter: acquiring 2 while having 1      
###                       (if >1 cooperative effect - if <1 competing)                      
### sigma21:              short-term interaction parameter: acquiring 1 while having 2      
###                       (if >1 cooperative effect - if <1 competing) 
### prop.immune:          proportion of immune cases (not used at the moment)
### n.seeds.1:            number of initial cases for path 1
### n.seeds.2:            number of initial cases for path 2
### rho[1]:               probability of being symptomatic for path 1
### rho[2]:               probability of being symptomatic for path 2
### alpha.as[1]:          relative infectiousness of asymptomatic cases (pathogen1)
### alpha.as[2]:          relative infectiousness of asymptomatic cases (pathogen2)
### netw:                 type of household network considered - Synthetic or ERGM
### n.vertex:             number of vertexes 
### n.networks:           number of simulated networks
### R.1:                  reproduction number path 1 (household R*)
### R.2:                  reproduction number path 2 (household R*)
### ratio.qhqg:           ratio transmission probability given household contact 
###                       over global contacts
### long.int[1]:          long-term interaction parameter: acquiring 2 while having  
###                       experienced (and recovered from) 1 
### long.int[2]:          long-term interaction parameter: acquiring 1 while having  
###                       experienced (and recovered from) 2
### pathogen[1]:          character variable identifying pathogen 1
### pathogen[2]:          character variable identifying pathogen 2
### contact.reduction:    parameter multiplying the household contact rate after home 
###                       isolation
### t.stop:               time at which simulations stop
### t.seed:               time of additional seeding
### behavior.change.1:    proportion of individuals changing behavior (home isolation) 
###                       after being infected with pathogen 1
### behavior.change.2:    proportion of individuals changing behavior (home isolation) 
###                       after being infected with pathogen 2
### reinfection:          boolean identifying whether someone can be re-infected with 
###                       the same pathogen (1 yes, 0 no)
### typeIC:               ID for different type of waning of immunity
### contact.reduction.TP: contact reduction value set to identify transmission rates 
###                       (household and global) linked to a specific R*
### behavior.change.1.TP: behavior change value (for pathogen 1) set to identify transmission 
###                       rates (household and global) linked to a specific R*
### behavior.change.2.TP: behavior change value (for pathogen 2) set to identify transmission 
###                       rates (household and global) linked to a specific R*
### het.vac:              boolean for heterologous effects (1 yes 0 no) - Not used currently
### t.imm.lim:            parameter to define the length of immunity that have the same overall 
###                       "effect" (area underneath the curve)
### decrease.gc:          decrease in the  number of global contact rates compared to baseline

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

sim.multipathogen <- function(HH.network, t2, t.seed){ #t2, lambda.g, sigma, prop.immune, n.seeds, rho, 
                              #inf.path.h, inf.path.g, alpha.as, long.int, pathogen, 
                              #contact.reduction,t.stop, t.seed, behavior.change, reinfection, 
                              #typeIC, het.vac, t.imm.lim){
#sim.multipathogen <- function(HH.network){

  #print("intitalize simulation")
  
  ######################################## STEP 1 ######################################## 
  # Create auxiliary objects 
  ########################################################################################
  
  netw.size <- network.size(HH.network)
  hh.id <- HH.network %v% "hh_id"
  hh.size <- HH.network %v% "hh_size"
  
  # status keeps track of the current status of rach individual. There is one table for each pathogen.
  status <- data.table(infected = rep(FALSE,netw.size),
                       recovered = rep(FALSE,netw.size),
                       time.infection = as.numeric(NA),
                       infected.by = as.numeric(NA),
                       severity = 0, # 1 symptomatic, 2 asymptomatic
                       time.symptom.onset = Inf,
                       immunity = 0, # 0 no immunity, 1 vaccinated
                       time.recovery = Inf)
  status <- list(status, status)
  
  # events contains the next event time for each possible event.
  events <- data.table(contact = Inf,
                       home.quarantine = as.numeric(Inf), # Wrap in list() to make it a list column
                       recovery = as.numeric(Inf),  # Wrap in list() to make it a list column
                       new.pathogen = t2,
                       new.seeding.1 = t.seed,
                       new.seeding.2 = t2 + t.seed)
  
  # table containing a logbook of all events during the simulation
  time.events <- setDT(data.frame(time = NA,
                           type = NA,
                           id = NA))  
  
  # transmission.parameters contains the transmission parameters for each individual for household
  # and global contacts, depending on their status.
  transmission.parameters <- setDT(data.frame(id = 1:netw.size, 
                                        q.1.h = as.numeric(0), 
                                        q.1.g = as.numeric(0), 
                                        q.2.h = as.numeric(0), 
                                        q.2.g = as.numeric(0), 
                                        contact.rate.h = sapply(1:netw.size, function(j) length(get.neighborhood(HH.network, j))),
                                        contact.rate.g = lambda.g, 
                                        susceptibility = 1))
  
  # home.quarantine contains information about each individual regarding home isolation. 
  home.quarantine <- setDT(data.frame(id = 1:netw.size,
                          quarantine = FALSE,
                          start = as.numeric(Inf),
                          end = as.numeric(Inf)))
  home.quarantine <- list(home.quarantine, home.quarantine)
  
  # next.contact contains the proposed time of the next contact for household and global contacts
  next.contact <- list(h = setDT(data.frame(id = 1:netw.size, 
                                     time.next.contact = Inf, 
                                     infectee.next.contact = NA)),
                      g = setDT(data.frame(id = 1:netw.size, 
                                     time.next.contact = Inf, 
                                     infectee.next.contact = NA)))
  
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
  
  current.time <- 0
  
  ######################################## STEP 2 ######################################## 
  # Introduce the first infections 
  ########################################################################################
  
  ### The individuals are randomly chosen in the population (among susceptibles).
  ### We assume that epidemics always start in household of size 2 
  ### (this is a random choice that decreases stochasticity).

  #print(paste0(current.time, ": Introducing the first pathogen"))
  
  first.cases.1 <- identify.initial.cases(pathogen[1], hh.size, status)
  for(j in first.cases.1){
    infection(infectee = j, path = 1, infector = 0, 
              current.time = current.time, rho = rho, behavior.change = behavior.change,
              env = environment(), inf.h = inf.h, inf.g = inf.g)
  }

  # variable to keep track of how often the acceptance rate was above 1
  err <- 0
  #Rt1 <- matrix(data = NA, nrow = 1, ncol = 2)
  #Rt2 <- matrix(data = NA, nrow = 1, ncol = 2)
  ### ?
  
  ##################################################################
  ### STEP 2
  ### Run the epidemic  
  ### Continue the epidemic as long as there are infected individuals or pathogen 2 still has to be introduced and as long
  ### as t.stop has not been reached.

  #while((sum(infectives)>0 & current.time<t.stop) | current.time<t2){ #while there are still infectives, we are within the t.stop
  while((sum(status[[1]]$infected + status[[2]]$infected) > 0 & current.time < t.stop) | current.time < t2){
    print(current.time)
    #print(paste0("Generate contacts and identify next events"))
    
    ### Phase 1: individuals that has to, propose a new social contact
    generate.next.contact(env = environment())
    
    ###  Phase 2: identify the next event: select the minimum time among the events that can occur
    next.event <- identify.next.event(env = environment())
    #print(next.contact$g, nrows = 101)
    
    #print(paste0("Next event: ", next.event))
    
    
    ### PHASE 3: the event happens
    if(next.event == "contact"){
      current.time <- events$contact
      #print(paste0("The next event is a contact at time ", round(current.time, 4)))
      ### STEP 1: select a contact and a contactee
      contacts <- list.next.contacts(env = environment())
      # if more than two contacts, sample one at random
      selected.contact <- contacts[sample(nrow(contacts),1),]
      
      ### STEP 2: "make" the contact
      if(selected.contact$type == "h"){ # for a household contact
        #print(paste0(selected.contact$id, " is making a household contact."))
        
        infectee.pool <- get.neighborhood(HH.network, selected.contact$id)
        
        # if multiple possibilities, choose one at random
        infectee <- ifelse(length(infectee.pool) != 0, 
                           infectee.pool[sample(length(infectee.pool),1)],
                           selected.contact$id)
        selected.contact$infectee.next.contact <- infectee
        next.contact$h$time.next.contact[selected.contact$id] <- Inf
      }else{ # for a global contact 
        infectee.pool <- setdiff(1:netw.size, c(selected.contact$id, get.neighborhood(HH.network, selected.contact$id)))
        infectee <- infectee.pool[sample(length(infectee.pool),1)]
        selected.contact$infectee.next.contact <- infectee
        next.contact$g$time.next.contact[selected.contact$id] <- Inf
      }
      
      #print(paste0("infectee pool: ", paste(infectee.pool, collapse = ", ")))
      #print(paste0("infectee: ", infectee))
      #print(paste0(round(current.time, 3), ": makes contact ", selected.contact$type))

      
      ### STEP 3: Infection
      # for each pathogen
      for(p1 in 1:length(pathogen)){
        # p2 is the indicator for the other pathogen
        p2 = setdiff(c(1:length(pathogen)), p1)
        # check if the infector is infected with disease p1
        if(status[[p1]]$infected[selected.contact$id] == TRUE){
          #print(paste0("Infector is infected with ", pathogen[p1]))
          # if infectee is not infected with p1, then compute short interaction
          if(status[[p1]]$infected[selected.contact$infectee.next.contact] == FALSE){
            #print(paste0("The infectee is not infected with pathogen ", pathogen[p1]))
            # if infectee is infected with the other disease, compute short term interaction
            short.inter <- ifelse(status[[p2]]$infected[selected.contact$infectee.next.contact] == TRUE,
                                  sigma[p1],
                                  1)
            # compute long-term interactions (i.e., the potential infectee already experienced the other infection, or that inefction)
            long.inter <- LLImmlev.basic(path.1 = p1,
                                         path.2 = p2,
                                         infectee = selected.contact$infectee.next.contact,
                                         t.imm.lim = t.imm.lim, env = environment())
            # select the transmission probability related to global or local contacts
            q <- ifelse(selected.contact$type == "g",
                        transmission.parameters[[paste0("q.", p1, ".g")]][selected.contact$id],
                        transmission.parameters[[paste0("q.", p1, ".h")]][selected.contact$id])
            # Acceptance rate is the probability that infection will follow from contact 
            # This is composed by q, short and long term interaction, and the infectiousness measure that describe how 
            # likely is that it will happen in that moment
            # if the acceptance rate is larger than 1, keep track of that. 
            acc.rate <- ifelse((selected.contact$type == "g" 
                                & home.quarantine[[p1]]$quarantine[selected.contact$infectee.next.contact] == TRUE) 
                               | status[[p1]]$infected[selected.contact$infectee.next.contact] == TRUE,
                               0,
                               InfMeasure(t = current.time - status[[p1]]$time.infection[selected.contact$id],
                                          path = p1) * short.inter * long.inter * q)
            #print(paste0("acceptance rate: ", acc.rate))
            if(acc.rate > 1){err <- err + 1}
            
            random <- runif(1)
            #print(paste0("random = ", random))
            if(random < acc.rate){
              #print("Infection!!!!")
              infection(infectee = selected.contact$infectee.next.contact,
                        path = p1, current.time = current.time, rho = rho,
                        infector = selected.contact$id, behavior.change = behavior.change,
                        env = environment(), inf.h = inf.h, inf.g = inf.g)
            }else{
              #print("Contact but no infection.")
              time.events <- rbind(time.events, list(current.time, paste0("unsuccesful contact type ", selected.contact$type), selected.contact$id))
            }
          }else{
            #print(paste0("infectee ", selected.contact$infectee.next.contact," is already infected with disease ", p1))
            #time.events <- rbind(time.events, list(current.time, paste0("unsuccesful contact type ", selected.contact$type), selected.contact$id))
          }
        }
      }
    }
    
    if(next.event=="home.quarantine"){ 
      current.time <- events$home.quarantine
      #print(paste0(current.time, ": home quarantine"))
      for(p in c(1:length(pathogen))){
        quarantines <- which(home.quarantine[[p]]$start == current.time)
        for(q in quarantines){
          home.quarantine[[p]][q, start := Inf]
          home.quarantine[[p]][q, stop := status[[p]]$time.recovery[q]]
          home.quarantine[[p]][q, quarantine := TRUE]
          
          next.contact$g[q, time.next.contact := Inf] 
          next.contact$h[q, time.next.contact := Inf] 
          transmission.parameters[q, contact.rate.h := transmission.parameters[q, contact.rate.h]*contact.reduction] 
        }
      }
    }
    
    if(next.event == "recovery"){ 
      # update the current time: 
      current.time <- events$recovery
      #print(paste0(current.time, ": recovery"))
      
      # for each disease
      for(p in c(1:2)){
        recoveries <- which(status[[p]]$time.recovery == events$recovery)
        for(r in recoveries){
          status[[p]][r, infected := FALSE]
          status[[p]][r, recovered := TRUE]
          status[[p]][r, time.recovery := Inf]
          time.events <- rbind(time.events, list(current.time, paste0("recovery from ", pathogen[p]), r))
          # come out of home quarantine
          home.quarantine[[p]][r, quarantine := FALSE]
          home.quarantine[[p]][r, start := Inf]
          home.quarantine[[p]][r, end := Inf]
          transmission.parameters$contact.rate.h[r] <- length(get.neighborhood(HH.network, r))
          # check if individual is not infected with other disease:
          if(status[[1]]$infected[r] + status[[2]]$infected[r] == 0){
            # remove all contacts that still had to happen in the future
            next.contact$g$time.next.contact[r] <- next.contact$h$time.next.contact[r] <- Inf
            # the person no longer needs to make a contact
          }
        }
      }
    }
    
    if(next.event == "new.pathogen"){ 
      # update the current time: 
      current.time <- events$new.pathogen
      #print(paste0(current.time, ": new pathogen"))
      
      # identify the individuals to be infected
      first.cases.2 <- identify.initial.cases(pathogen[2], hh.size, status)
      # infect those individuals
      for(j in first.cases.2){
        infection(infectee = j, path = 2, infector = 0, current.time = current.time, 
                  rho = rho, inf.h = inf.h, inf.g = inf.g, behavior.change = behavior.change, 
                  env = environment())
      }
      # the event has been introduced, the time needs to be set to Inf. 
      events[, new.pathogen := Inf]
    }
    
    if(next.event == "new.seeding.1"){
      # update the current time: 
      current.time <- events$new.seeding.1
      #print(paste0(current.time, ": new seeding 1"))
      
      # identify the susceptible individuals
      not.infected <- which(status[[1]]$infected!=TRUE)
      
      if (n.seeds[1] < length(not.infected)){ # if there are more susceptibles than seeds
        new.cases <- sample(not.infected, n.seeds[1])
        for(j in new.cases){
          infection(infectee = j, path = pathogen[1], infector = 0,
                    current.time = current.time, rho = rho, inf.h = inf.h,
                    inf.g = inf.g, behavior.change = behavior.change, env = environment())
        }
      }
      # update the event times: 
      #events$new.seeding.1 <- current.time+t.seed
      events$new.seeding.1 <- Inf
    }
    
    if(next.event == "new.seeding.2"){
      # update the current time: 
      current.time <- events$new.seeding.2
      #print(paste0(current.time, ": new seeding 2"))
      
      # identify the susceptible individuals
      not.infected <- which(status[[2]]$infected!=TRUE)
      
      if (n.seeds.2 < length(not.infected)){ # if there are more susceptibles than seeds
        new.cases <- sample(not.infected, n.seeds[2])
        for(j in new.cases){
          infection(infectee = j, path = pathogen[1], infector = 0,
                    current.time = current.time, rho = rho, inf.h = inf.h,
                    inf.g = inf.g, behavior.change = behavior.change, env = environment())
        }
      }
      # update the event times: 
      events$new.seeding.2 <<- Inf
    }
  }
  
  time.events <- na.omit(time.events) # delete the first empty row from data frame
  
  # compute some summary measures that will be given as output
  #C1 <- n.seeds.1
  #Y1 <- n.seeds.1
  #if(t2 > 0){
  #  C2 <- 0
  #  Y2 <- 0
  #}else{
  #  C2 <- n.seeds.2
  #  Y2 <- n.seeds.2
  #}
  #last.day <- round(max(time.events$time))
  
  #for (i in 1:last.day){
  #  temp.time <- setdiff(which(time.events[,1] > i), which(time.events[,1] > i+1))
  #  temp.inf.1<-c(which(time.events[temp.time,2]==1.1),which(time.events[temp.time,2]==1.2))
  #  temp.inf.2<-c(which(time.events[temp.time,2]==2.1),which(time.events[temp.time,2]==2.2))
  #  temp.time.1<-setdiff(1:length(time.events[,1]),which(time.events[,1]>i+1))
  #  C1<- c(C1,length((which(time.events[temp.time.1,2]==1.1)))+length((which(time.events[temp.time.1,2]==1.2)))-length((which(time.events[temp.time.1,2]==-1))))
  #  C2<- c(C2,length((which(time.events[temp.time.1,2]==2.1)))+length((which(time.events[temp.time.1,2]==2.2)))-length((which(time.events[temp.time.1,2]==-2))))
  #  Y1<- c(Y1,length(temp.inf.1))
  #  Y2<-c(Y2,length(temp.inf.2))
  #}
  
  #Fs1<-length(which(time.events[,2]==1.1))+length(which(time.events[,2]==1.2))
  #Fs2<-length(which(time.events[,2]==2.1))+length(which(time.events[,2]==2.2))
  
  #epi.details<-data.frame("Days"=0:last.day, "Incidence1"=Y1,"Incidence2"=Y2, "Prevalence1"=C1,"Prevalence2"=C2)
  #FinalSize<-data.frame("FinalSize1"=Fs1,"FinalSize2"=Fs2)
  #PeakIncidence<-data.frame("PeakIncidence1"=max(epi.details$Incidence1),"TimePeakIncidence1"=which(epi.details$Incidence1==max(epi.details$Incidence1))[1],"PeakIncidence2"=max(epi.details$Incidence2),"TimePeakIncidence2"=which(epi.details$Incidence2==max(epi.details$Incidence2))[1] )
  #PeakPrevalence<-data.frame("PeakPrevalence1"=max(epi.details$Prevalence1),"TimePeakPrevalence1"=which(epi.details$Prevalence1==max(epi.details$Prevalence1))[1],"PeakPrevalence2"=max(epi.details$Prevalence2),"TimePeakPrevalence2"=which(epi.details$Prevalence2==max(epi.details$Prevalence2))[1] )
  #return(list(time.events=time.events, status.matrix.1=status.matrix.1, status.matrix.2=status.matrix.2,epi.details=epi.details, FinalSize=FinalSize, PeakIncidence=PeakIncidence, PeakPrevalence=PeakPrevalence,Rt1=Rt1,Rt2=Rt2))
}

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

InfMeasure <- function(t, path){
  inf.measure = 0
  if(pathogen[path] %in% c("COVID-19", "DELTA", "OMICRON")){
    inf.measure = dgamma(t, shape = 12, rate = 2.08)/(pgamma(12, shape = 12, rate = 2.08))
  }
  if(pathogen[path] == "FLU-A"){
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
    inf.measure = 1.001592*dgamma(t, shape = 4.604016, scale = 0.5922317)
    #return(dgamma(t,shape = 3.5, rate = 1.15)/(pgamma(6.24,shape = 3.5,rate = 1.15)))
  }
  if(pathogen[path] == "RSV"){
    inf.measure = dgamma(t, shape = 15, rate = 2.6)/(pgamma(12, shape = 15, rate = 2.6))
  }
  return(inf.measure)
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

LLImmlev.basic <- function(path.1,
                           path.2,
                           infectee, 
                           t.imm.lim, env){ #pathogen.v1 is the infection the infectee might catch
  
  value <- 1
  if(env$status[[path.2]]$recovered[infectee] == TRUE){
    # how long ago did the individual recover from the other disease?
    t.since.infection <- env$current.time - (env$status[[path.2]]$time.infection[infectee] 
                                         + infectious.period.length(pathogen[path.2]))
    if(t.since.infection < t.imm.lim){ # if the infection was recent enough
      if(typeIC == 1){
        value <- long.int[path.1]
      }
      if(typeIC == 2){
        value <- (t.since.infection/t.imm.lim)
      }
      if(typeIC == 3){
        value <- 0.25 + (t.since.infection/(2*t.imm.lim))
      }
      if(typeIC == 4){
        value <- 2.5 - t.since.infection/10
      }
      if(typeIC == 5){
        value <- long.int[path.1]
      }
      if(typeIC == 6){
        value <- 3 - t.since.infection/5
      }
    }  
  }
  return(value)
}


################################################
### This function computes RT                ###
### started but haven't look at this anymore ###
################################################

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

##########################################################################
### This function choses the initial infections for a certain pathogen ###
##########################################################################

identify.initial.cases = function(path, hh.size, status){
  p <- ifelse(path == pathogen[1], 1, 2)
  # potential seeds are all individuals from households of size 2
  potential.seeds <- which(hh.size == 2)
  # from those potential seeds that are not yet infected with path, n.seeds individuals are randomly sampled
  first.cases <- sample(potential.seeds[which(status[[p]][potential.seeds, "infected"] == FALSE)], n.seeds[p])
  return(first.cases)
}

#########################################################################
### This function does what has to happen during an event "infection" ###
#########################################################################

infection <- function(infectee, path, infector, current.time, rho, inf.h, inf.g, behavior.change, env){
  #print("Start infection now.")
  # index for the pathogen
  p <- path
  #print(pathogen[p])
  # create a deep copy to prevent both data tables in the list to be adapted.
  status_copy <- copy(env$status[[p]]) #create a deep copy to prevent both data tables in the list to be adapted.
  home.quarantine_copy <- copy(env$home.quarantine[[p]])
  
  status_copy[infectee, infected := TRUE]
  status_copy[infectee, time.infection := current.time]
  status_copy[infectee, time.recovery := current.time + infectious.period.length(pathogen = pathogen[p])]
  status_copy[infectee, infected.by := infector]
  #print(status_copy[infectee])
  
  if(runif(1) < rho[p]){ # if the individual is symptomatic
    #print("symptomatic infection")
    env$transmission.parameters[infectee, (paste0("q.", p, ".h")) := inf.h[p]]
    env$transmission.parameters[infectee, (paste0("q.", p, ".g")) := inf.g[p]]
    #print(transmission.parameters[infectee])
    status_copy[infectee, severity := 1]
    status_copy[infectee, time.symptom.onset := current.time+incubation.period(pathogen[p])]
    #print(status_copy[infectee])
    if(runif(1) < behavior.change[p]){
      #print("behavior change")
      home.quarantine_copy[infectee, quarantine := TRUE]
      home.quarantine_copy[infectee, start := status_copy[infectee, time.symptom.onset]]
      home.quarantine_copy[infectee, end := status_copy[infectee, time.recovery]]
      #print(home.quarantine_copy[infectee])
    }
    env$time.events <- rbind(env$time.events, list(current.time, paste0("symptomatic infection with ", pathogen[p]), infectee))
    #print(time.events)
  }else{
    #print("asymptomatic infection")
    
    env$transmission.parameters[infectee, (paste0("q.", p, ".h")) := inf.h[p] * alpha.as[p]]
    env$transmission.parameters[infectee, (paste0("q.", p, ".g")) := inf.g[p] * alpha.as[p]]
    #print(transmission.parameters[infectee])
    
    status_copy[infectee, severity := 2]
    #print(status_copy[infectee])
    
    env$time.events <- rbind(env$time.events, list(current.time, paste0("asymptomatic infection with ", pathogen[p]), infectee))
    #print(time.events)
  }
  env$status[[p]] <- status_copy
  env$home.quarantine[[p]] <- home.quarantine_copy
}

#################################################
### This function generates the next contacts ###
#################################################

generate.next.contact = function(env){
  infected <- which(env$status[[1]]$infected == TRUE | env$status[[2]]$infected == TRUE)
  household.contacts <- setdiff(infected, which(env$next.contact$h$time.next.contact < Inf))
  global.contacts <- setdiff(setdiff(infected, which(env$next.contact$g$time.next.contact < Inf)),
                            c(env$home.quarantine[[1]][env$home.quarantine[[1]]$quarantine==TRUE,]$id,
                                  env$home.quarantine[[2]][env$home.quarantine[[2]]$quarantine==TRUE,]$id))
                   
  for(j in household.contacts){
    env$next.contact$h[j,time.next.contact := ifelse(env$transmission.parameters[j, contact.rate.h] != 0, # check whether there are any household members
                                                 rexp(1, env$transmission.parameters[j, contact.rate.h]) + env$current.time,
                                                 Inf)]
  }
  
  for(j in global.contacts){
    env$next.contact$g[j, time.next.contact := rexp(1, env$transmission.parameters$contact.rate.g[j]) + env$current.time]
  }
}

###############################################
### This function identifies the next event ###
###############################################

identify.next.event = function(env){
  # identify the next contact among all household and global contacts
  temp1 <- min(c(env$next.contact$h[, time.next.contact], env$next.contact$g[, time.next.contact]), na.rm = TRUE)
  env$events[, contact := temp1]
  
  # identify the next home.quarantine time for both pathogens
  temp2 <- min(c(env$home.quarantine[[1]]$start, env$home.quarantine[[2]]$start), na.rm = TRUE)
  env$events[, home.quarantine := temp2]
 
  # identify the next recovery time for both pathogens
  temp3 <- min(c(env$status[[1]]$time.recovery, env$status[[2]]$time.recovery), na.rm = TRUE)
  env$events[, recovery := temp3]
  
  # select the next event(s). If more than one event: sample one at random
  next.event <- sample(names(env$events)[which.min(env$events)], 1)
  #if(length(next.event) > 1){
  #  next.event <- sample(next.events ,1)
  #}
  return(next.event)
}

###################################################################################
### This function returns a data set with all contacts that need to happen next ###
###################################################################################

list.next.contacts = function(env){
  contacts.g <- env$next.contact$g[env$next.contact$g$time.next.contact == (min(env$next.contact$g$time.next.contact, 
                                                                        env$next.contact$h$time.next.contact)),]
  try(contacts.g$type <- "g", silent = TRUE)
  contacts.h <- env$next.contact$h[env$next.contact$h$time.next.contact == (min(env$next.contact$g$time.next.contact, 
                                                                        env$next.contact$h$time.next.contact)),]
  try(contacts.h$type <- "h", silent = TRUE)
  contacts <- rbind(contacts.g, contacts.h)
  return(contacts)
}















