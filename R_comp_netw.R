# Function to compute the value of the reproduction number in a 2-level household mixing. 
#Constant infectious period

R0.computation.StaticNetw<-function(HH.network,beta.g,nSim, beta.h,prob.asym,asymp.rel.inf,lambda.h){
  mu<-1
  hh.id<-HH.network %v% "hh_id"
  hh.size<-HH.network %v% "hh_size"
  n<-length(hh.id)
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
    n.asympt<-0
    q.h<-beta.h/lambda.h
      for (w in unique(hh.id)){
        #print(w)
        hh.data<-data.frame("members"= which(hh.id==w),"id"=1:length(which(hh.id==w)),"status"=0,"recovery"=Inf, "index.contact"=0, "betah"=0)
        #index case
        for (r in hh.data$id){
          ifelse(length(get.neighborhood(HH.network, hh.data$members[r]))>0,hh.data$betah[r]<-(length(get.neighborhood(HH.network, hh.data$members[r]))*q.h),hh.data$betah[r]<-1/rexp(1,1/exp(100))) 
          if (runif(1)<prob.asym){
            hh.data$betah[r]<-hh.data$betah[r]*asymp.rel.inf
            n.asympt<-n.asympt+1
          }
        }
        primary<-sample(1:length(hh.data$members),1)
        hh.data$status[primary]<-1 
        hh.data$recovery[primary]<-mu
        hh.data$index.contact[primary]<-1
        contact.time<-data.frame("id"=hh.data$members,"pr.ctc"=rep(NA,length(hh.data$members)),"pr.infectee"=rep(NA,length(hh.data$members)))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
        current.time<-0
        while ((length(which(hh.data$status==1))>0)) { # till there is at least one infectious individual
          for (i in which(hh.data$index.contact==1) ){ # for all the individuals that has to propose a new contact
            temp.contact.time<-rexp(1,hh.data$betah[i])+current.time
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
    
    if (n.asympt==n){
      R0[j]<-(beta.g*asymp.rel.inf)*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
    }
    if (n.asympt==0){
      R0[j]<-(beta.g)*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
    }
    if (n.asympt!=n & n.asympt!=0){
      R0[j]<-((beta.g*n/(n- n.asympt))+ (beta.g*asymp.rel.inf*n.asympt/n))*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
    }  

    
    print(j)
  }
  return(R0)
}

R0.comp<-function(ratio_hhgl,tol,R.rif,HH.network,nSim,prob.asym,asymp.rel.inf,lambda.h){
  beta.g<-1
  beta.h<-ratio_hhgl*beta.g
  beta.g.tempm<-0
  beta.g.tempM<-10
  R.temp<-NULL
  for (i in 1:nSim){
    temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
    R.temp<-c(R.temp,R0.computation.StaticNetw(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1,prob.asym=prob.asym,asymp.rel.inf=asymp.rel.inf,lambda.h = lambda.h))
  }
  R.temp<-mean(R.temp)
  while (abs(mean(R.temp)-R.rif)>tol){
    R.temp<-NULL
    for (i in 1:nSim){
      temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
      R.temp<-c(R.temp,R0.computation.StaticNetw(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1,prob.asym=prob.asym,asymp.rel.inf=asymp.rel.inf,lambda.h = lambda.h))
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

R0.computation.Inf<-function(HH.network,beta.g,nSim, beta.h,prob.asym,asymp.rel.inf,lambda.h,pathogen,ctc.dec,compl){
  mu<-infectious.period.length(pathogen = pathogen)
  meanIP<-mean.ip(pathogen = pathogen)
  hh.id<-HH.network %v% "hh_id"
  hh.size<-HH.network %v% "hh_size"
  n<-length(hh.id)
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
    n.asympt<-0
    q.h<-beta.h/lambda.h
    for (w in unique(hh.id)){
      #print(w)
      hh.data<-data.frame("members"= which(hh.id==w),"id"=1:length(which(hh.id==w)),"status"=0,"recovery"=Inf, "index.contact"=0, "betah"=0, "SO"=Inf, "ToI"=Inf)
      #index case
      for (r in hh.data$id){
        ifelse(length(get.neighborhood(HH.network, hh.data$members[r]))>0,hh.data$betah[r]<-(length(get.neighborhood(HH.network, hh.data$members[r]))*q.h),hh.data$betah[r]<-1/rexp(1,1/exp(100))) 
        if (runif(1)<prob.asym){
          hh.data$betah[r]<-hh.data$betah[r]*asymp.rel.inf
          n.asympt<-n.asympt+1
        }
      }
      primary<-sample(1:length(hh.data$members),1)
      hh.data$status[primary]<-1 
      hh.data$recovery[primary]<-mu
      hh.data$index.contact[primary]<-1
      hh.data$SO[primary]<-incubation.period(pathogen = pathogen)
      hh.data$ToI[primary]<-0
      contact.time<-data.frame("id"=hh.data$members,"pr.ctc"=rep(NA,length(hh.data$members)),"pr.infectee"=rep(NA,length(hh.data$members)))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
      current.time<-0
      events<-data.frame(NextCtc        = Inf,
                         SymptOns       = Inf,
                         Recovery       = Inf)
      while ((length(which(hh.data$status==1))>0)) { # till there is at least one infectious individual
        for (i in which(hh.data$index.contact==1) ){ # for all the individuals that has to propose a new contact
          temp.contact.time<-rexp(1,hh.data$betah[i])+current.time
           hh.members.contacted<-get.neighborhood(HH.network,hh.data$members[i], type = "out")
          hh.data$index.contact[i]<-0
          ifelse(length(hh.members.contacted)>0,contact.time$pr.infectee[i] <-sample(hh.members.contacted,1),contact.time$pr.infectee[i] <-NA)
          if (length(hh.members.contacted)==1){contact.time$pr.infectee[i]<-hh.members.contacted}
          ifelse(length(hh.members.contacted)>0,contact.time$pr.ctc[i]<-temp.contact.time ,contact.time$pr.ctc[i]<-NA)
        }
        #computation of the next event
        ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,events$NextCtc<-min(contact.time$pr.ctc, na.rm = T),events$NextCtc<-Inf)
        events$Recovery<-min(hh.data$recovery, na.rm = T)
        events$SymptOns<-min(hh.data$SO, na.rm = T)
        
        next.evts<-colnames(events)[min(events)==events]
        if (length(next.evts)>1){
          next.evts<-sample(next.evts,1)
        }
        #next event is an infection
        
        if (next.evts=="NextCtc"){
          current.time<-events$NextCtc
          infector<-which(contact.time$pr.ctc ==current.time)
          infectee<-hh.data$id[which(hh.data$members==contact.time$pr.infectee[infector])]
          if (hh.data$status[infectee]==0 & runif(1)<InfMeasure(t=current.time-hh.data$ToI[infector], pathogen = pathogen)){
            hh.data$status[infectee]<-1
            hh.data$recovery[infectee]<-current.time+mu
            hh.data$SO[infectee]<-current.time+incubation.period(pathogen = pathogen)
            hh.data$ToI[infectee]<-current.time
            hh.data$index.contact[infectee]<-1
            hh.data$index.contact[infector]<-1
            contact.time[infector,2:3]<-NA
          }else{
            hh.data$index.contact[infector]<-1
            contact.time[infector,2:3]<-NA
          }
          
        }
        if (next.evts=="Recovery"){
          current.time<-events$Recovery
          recovered<-which(hh.data$recovery==current.time)
          hh.data$recovery[recovered]<-Inf
          hh.data$status[recovered]<--1
          contact.time[recovered,2:3]<-rep(NA,2)
        }
        if (next.evts=="SymptOns"){
          current.time<-events$SymptOns
          symptomatic<-which(hh.data$SO==current.time)
          hh.data$SO[symptomatic]<-Inf
          hh.data$betah[symptomatic]<-hh.data$betah[symptomatic]*ctc.dec
          contact.time[symptomatic,2:3]<-rep(NA,2)
          hh.data$index.contact[symptomatic]<-1
        }
  
      }
      AR[[hh.size[hh.data$members[1]]]]<-c(AR[[hh.size[hh.data$members[1]]]],length(which(hh.data$status==-1)))
      # sar - NA when no symptomatic infection are register
    }
    
    ar<-0
    for (s in 1:max(unique(hh.size))){
      ar[s]<-ifelse(length(AR[[s]])>1,mean(AR[[s]][-1]),0) 
    }
    
    if (compl==0){
      if (n.asympt==n){
        R0[j]<-(beta.g*mu*asymp.rel.inf)*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
      }
      if (n.asympt==0){
        R0[j]<-(beta.g*mu)*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
      }
      if (n.asympt!=n & n.asympt!=0){
        R0[j]<-((beta.g*mu*(n- n.asympt)/n)+ (beta.g*mu*asymp.rel.inf*n.asympt/n))*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
      }  
      
    }else{
      if (n.asympt==n){
        R0[j]<-(beta.g*asymp.rel.inf)*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
      }
      if (n.asympt==0){
        R0[j]<-(beta.g*meanIP*compl+beta.g*mu*(1-compl)) *(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
      }
      if (n.asympt!=n & n.asympt!=0){
        R0[j]<-(((beta.g*meanIP*compl+beta.g*mu*(1-compl))*(n- n.asympt)/n)+ (beta.g*mu*asymp.rel.inf*n.asympt/n))*(sum(ar*(h.n)*(1:max(unique(hh.size)))))/mu.h
      }  
    }
    print(j)
  }
  return(R0)
}

R0.comp.Inf<-function(ratio_hhgl,tol,R.rif,HH.network,nSim,prob.asym,asymp.rel.inf,lambda.h,pathogen,ctc.dec,compl){
  mu<-infectious.period.length(pathogen = pathogen)
  beta.g<-1
  beta.h<-ratio_hhgl*beta.g
  beta.g.tempm<-0
  beta.g.tempM<-10
  R.temp<-NULL
  for (i in 1:nSim){
    temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
    R.temp<-c(R.temp,R0.computation.Inf(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1,prob.asym=prob.asym,asymp.rel.inf=asymp.rel.inf,lambda.h = lambda.h,pathogen = pathogen,ctc.dec = ctc.dec, compl = compl))
  }
  R.temp<-mean(R.temp)
  while (abs(mean(R.temp)-R.rif)>tol){
    R.temp<-NULL
    for (i in 1:nSim){
      temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
      R.temp<-c(R.temp,R0.computation.Inf(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1,prob.asym=prob.asym,asymp.rel.inf=asymp.rel.inf,lambda.h = lambda.h,pathogen = pathogen,ctc.dec = ctc.dec, compl = compl))
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
  transm.prms<-data.frame("beta.g"=beta.g*mu, "beta.h"=beta.h)
  return(transm.prms)
}

Comp.R0<-function(HH.network,nSim,transm.prms,prob.asym,asymp.rel.in,pathogen,ctc.dec,compl){
  beta.g<-transm.prms$beta.g
  beta.h<-transm.prms$beta.h
  R.temp<-NULL
  for (i in 1:nSim){
    temp.HH.netw<-HH.network[[sample(1:length(HH.network),1)]]
    R.temp<-c(R.temp,R0.computation.Inf(HH.network = temp.HH.netw, beta.g = beta.g, beta.h = beta.h, nSim = 1,prob.asym=prob.asym,asymp.rel.inf=asymp.rel.inf,lambda.h = lambda.h,pathogen = pathogen,ctc.dec = ctc.dec, compl = compl))
  }
  return(R.temp)
}

incubation.period<-function(pathogen){
  if (pathogen=="COVID-19" | pathogen == "DELTA" | pathogen== "OMICRON"){
    return(rlnorm(1,meanlog = log(5.2), sdlog = log(1.7)))    
  }
  if (pathogen=="FLU-A"){
    return(rlnorm(1,meanlog = log(1.4), sdlog = log(1.51))) 
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

mean.ip<-function(pathogen){
  if (pathogen=="COVID-19" | pathogen == "DELTA" | pathogen== "OMICRON"){
    return(5.288462)    
  }
  if (pathogen=="FLU-A"){
    return(2) 
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
    return(15)
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
    return(dgamma(t,shape = 12, rate = 2.08)/ (pgamma(15,shape = 12,rate = 2.08)))
  }
  if (pathogen=="FLU-A"){
    # Setting infectiousness measure according to Carrat et al. (2008) for H1N1
    #VL<-data.frame(x=0:8,y=c(0,1.75,3,2.5,1.8,1.25,0.75,0.5,0)) 
    #vl.flu<-nlsLM(y~a*dgamma(x=x,shape = s1,scale = sc1),start = list(a=10,s1=1.5,sc1=1.5),data = VL, weights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1))
    #prms<-vl.flu$m$getPars()
    #print(prms)
    #plot(seq(0,8,0.1),prms[1]*dgamma(x=seq(0,8,0.1),shape = prms[2],scale = prms[3]),ylim = c(0,3.1))
    #points(VL, col="red")
    #f.vl<-function(t){
    #  return(prms[1]*dgamma(t,shape = prms[2],scale = prms[3]))
    #}
    #k<-integrate(f.vl,lower = 0,upper = 8)
    #prms[1]/ k$value
    return(1.014939*dgamma(t,shape = 3.2930139,scale = 0.9533298))
    
    
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
