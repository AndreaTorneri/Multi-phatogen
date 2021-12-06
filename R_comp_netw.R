# Function to compute the value of the reproduction number in a 2-level household mixing. 
#Constant infectious period





R0.computation.StaticNetw<-function(HH.network,beta.g,nSim, beta.h){
  mu<-1
  hh.id<-HH.network %v% "hh_id"
  hh.size<-HH.network %v% "hh_size"
  m<- length(unique(hh.id)) # number of households
  m.n<-rep(0,max(unique(hh.size))) # number of household of a specific size (position 1 size 1, position 2 size 2, ....)
  R0<-0
  for (i in 1:length(m.n)){
  m.n[i]<-length(unique(hh.id[which(hh.size==i)]))
}  
  h.n<-m.n/m #proportion of household of a specific size
  mu.h<-sum(h.n*(1:length(h.n)))
    
  #divide network in sizes
  for (j in 1:nSim){
    AR<-list()
    for (s in 1:max(unique(hh.size))) {
      AR[[s]]<-0  
    }
      for (w in unique(hh.id)){
        #print(w)
        hh.data<-data.frame("members"= which(hh.id==w),"id"=1:length(which(hh.id==w)),"status"=0,"recovery"=Inf, "index.contact"=0)
        #index case
        primary<-sample(1:length(hh.data$members),1)
        hh.data$status[primary]<-1 
        hh.data$recovery[primary]<-mu
        hh.data$index.contact[primary]<-1
        contact.time<-data.frame("id"=hh.data$members,"pr.ctc"=rep(NA,length(hh.data$members)),"pr.infectee"=rep(NA,length(hh.data$members)))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
        current.time<-0
        while ((length(which(hh.data$status==1))>0)) { # till there is at least one infectious individual
          for (i in which(hh.data$index.contact==1) ){ # for all the individuals that has to propose a new contact
            temp.contact.time<-rexp(1,beta.h)+current.time
            hh.members.contacted<-get.neighborhood(HH.network,hh.data$members[i], type = "out")
            hh.data$index.contact[i]<-0
              ifelse(length(hh.members.contacted)>0,contact.time$pr.infectee[i] <-sample(hh.members.contacted,1),contact.time$pr.infectee[i] <-NA)
              if (length(hh.members.contacted)==1){contact.time$pr.infectee[i]<-hh.members.contacted}
              ifelse(length(hh.members.contacted)>0,contact.time$pr.ctc[i]<-temp.contact.time ,contact.time$pr.ctc[i]<-NA)
            }
          #computation of the next event
          ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_h<-min(contact.time$pr.ctc, na.rm = T),T_h<-Inf)
          R_a<-min(hh.data$recovery, na.rm = T)
          
          #next event is an infection
          if (T_h<R_a){
            # a household contact
              current.time<-T_h
              infector<-which(contact.time$pr.ctc ==T_h)
              infectee<-hh.data$id[which(hh.data$members==contact.time$pr.infectee[infector])]
              if (hh.data$status[infectee]==0 ){
                hh.data$status[infectee]<-1
                hh.data$recovery[infectee]<-current.time+mu
                hh.data$index.contact[infectee]<-1
                hh.data$index.contact[infector]<-1
                contact.time$pr.ctc[infector]<-NA
              }else{
              hh.data$index.contact[infector]<-1
              contact.time[infector,2:3]<-NA
            }
            #Phase 2.3 a recovery occurs  
          }else{
            current.time<-R_a
            recovered<-which(hh.data$recovery==R_a)
            hh.data$recovery[recovered]<-Inf
            hh.data$status[recovered]<--1
            contact.time[recovered,2:3]<-rep(NA,2)
          }
          }
        AR[[hh.size[hh.data$members[1]]]]<-c(AR[[hh.size[hh.data$members[1]]]],length(which(hh.data$status==-1)))
        # sar - NA when no symptomatic infection are register
        }
    
    ar<-0
    for (s in 1:max(unique(hh.size))){
    ar[s]<-ifelse(length(AR[[s]])>1,mean(AR[[s]][-1]),0) 
    }
    
  R0[j]<-beta.g*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
    print(j)
  }
  return(R0)
}

R0.comp<-function(ratio_hhgl,tol,R.rif,HH.network,nSim){
  beta.g<-1
  beta.h<-ratio_hhgl*beta.g
  beta.g.tempm<-0
  beta.g.tempM<-10
  R.temp<-NULL
  for (i in 1:nSim){
    temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
    R.temp<-c(R.temp,R0.computation.StaticNetw(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1))
  }
  R.temp<-mean(R.temp)
  while (abs(mean(R.temp)-R.rif)>tol){
    R.temp<-NULL
    for (i in 1:nSim){
      temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
      R.temp<-c(R.temp,R0.computation.StaticNetw(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1))
    }
    R.temp<-mean(R.temp)
    if (mean(R.temp)>R.rif){
      beta.g.tempM<-beta.g
      beta.g<-runif(1,min = beta.g.tempm,max = beta.g)
      beta.h<-beta.g*ratio_hhgl
    }else{
      beta.g.tempm<-beta.g
      beta.g<-runif(1,min = beta.g, max = beta.g.tempM)
      beta.h<-beta.g*ratio_hhgl
    }
    print(c(R.rif,mean(R.temp),abs(mean(R.temp)-R.rif)))
  }
  transm.prms<-data.frame("beta.g"=beta.g, "beta.h"=beta.h)
  return(transm.prms)
}

R0.comp.asy<-function(ratio_hhgl,tol,R.rif,HH.network,nSim,IPL){
  beta.g<-1
  beta.h<-ratio_hhgl*beta.g
  beta.g.tempm<-0
  beta.g.tempM<-10
  R.temp<-R0.computation.StaticNetw(HH.network = HH.network, beta.g = beta.g, beta.h = beta.h, nSim = nSim)
  while (abs(mean(R.temp)-R.rif)>tol){
    R.temp<-R0.computation.StaticNetw(HH.network = HH.network, beta.g = beta.g, beta.h = beta.h, nSim = nSim)
    if (mean(R.temp)>R.rif){
      beta.g.tempM<-beta.g
      beta.g<-runif(1,min = beta.g.tempm,max = beta.g)
      beta.h<-beta.g*ratio_hhgl
    }else{
      beta.g.tempm<-beta.g
      beta.g<-runif(1,min = beta.g, max = beta.g.tempM)
      beta.h<-beta.g*ratio_hhgl
    }
    print(c(R.rif,mean(R.temp),abs(mean(R.temp)-R.rif)))
  }
  transm.prms<-data.frame("beta.g"=beta.g, "beta.h"=beta.h)
  return(transm.prms)
}


avg.lambda.within<-function(HH.networks){
  lamb.with<-NULL
  for (i in 1:length(HH.networks)){
    temp.net<-HH.networks[[i]]
    vn<-temp.net  %v% "vertex.names"
    temp.lh<-NULL
    for (j in 1:length(vn)){
      temp.lh<-c(temp.lh, length(get.neighborhood(temp.net,j)))
    }
    lamb.with<-c(lamb.with,mean(temp.lh))
  }
  return(mean(lamb.with))
}