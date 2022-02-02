
creat.synth.netw<-function(n.vertex,n.networks,density_by_hh_size){
  # Computing the ratio q_h/q_g
  library(network)
  
  HH.netw<-list()
  n<-n.vertex
  # Compute the synthetic scenarios
  max.hhsz<-7
  data.hhnetw<-density_by_hh_size[1:max.hhsz,]
  prob.hhsize<-data.hhnetw$household_size_freq/sum(data.hhnetw$household_size_freq)
  count<-0
  hh.size.synth<-data.frame("Size"=1:max.hhsz,"Freq"=0)
  
  for (w in 1:n.networks){
    while(count<n){
      k<-sample(1:max.hhsz,1,prob = prob.hhsize)
      hh.size.synth$Freq[which(hh.size.synth$Size==k)]<-hh.size.synth$Freq[which(hh.size.synth$Size==k)]+1
      count<-count+k
    }
    
    adj.mat<-matrix(data=0,nrow = count, ncol = count)
    hh.size.vect<-rep(1,count)
    temp21<-(hh.size.synth$Freq[1]+1)
    temp22<-hh.size.synth$Freq[1]+2*hh.size.synth$Freq[2]
    for(j in seq(temp21,temp22,by=2) ){
      adj.mat[j,j+1]<-1
      adj.mat[j+1,j]<-1
      hh.size.vect[j:(j+1)]<-2
    }
    
    temp31<-temp22+1
    temp32<-hh.size.synth$Freq[1]+2*hh.size.synth$Freq[2]+3*hh.size.synth$Freq[3]
    for(j in seq(temp31,temp32,by=3) ){
      adj.mat[j,j+1]<-1
      adj.mat[j,j+2]<-1
      adj.mat[j+1,j]<-1
      adj.mat[j+1,j+2]<-1
      adj.mat[j+2,j]<-1
      adj.mat[j+2,j+1]<-1
      hh.size.vect[seq(j,j+2,1)]<-3
    }
    
    temp41<-temp32+1
    temp42<-hh.size.synth$Freq[1]+2*hh.size.synth$Freq[2]+3*hh.size.synth$Freq[3]+4*hh.size.synth$Freq[4]
    for(j in seq(temp41,temp42,by=4) ){
      adj.mat[j,j+1]<-1
      adj.mat[j,j+2]<-1
      adj.mat[j,j+3]<-1
      adj.mat[j+1,j]<-1
      adj.mat[j+1,j+2]<-1
      adj.mat[j+1,j+3]<-1
      adj.mat[j+2,j]<-1
      adj.mat[j+2,j+1]<-1
      adj.mat[j+2,j+3]<-1
      adj.mat[j+3,j+1]<-1
      adj.mat[j+3,j+2]<-1
      adj.mat[j+3,j]<-1
      hh.size.vect[j:(j+3)]<-4
    }
    
    temp51<-temp42+1
    temp52<-hh.size.synth$Freq[1]+2*hh.size.synth$Freq[2]+3*hh.size.synth$Freq[3]+4*hh.size.synth$Freq[4]+5*hh.size.synth$Freq[5]
    for(j in seq(temp51,temp52,by=5) ){
      adj.mat[j,j+1]<-1
      adj.mat[j,j+2]<-1
      adj.mat[j,j+3]<-1
      adj.mat[j,j+4]<-1
      adj.mat[j+1,j]<-1
      adj.mat[j+1,j+2]<-1
      adj.mat[j+1,j+3]<-1
      adj.mat[j+1,j+4]<-1
      adj.mat[j+2,j+1]<-1
      adj.mat[j+2,j]<-1
      adj.mat[j+2,j+3]<-1
      adj.mat[j+2,j+4]<-1
      adj.mat[j+3,j+1]<-1
      adj.mat[j+3,j+2]<-1
      adj.mat[j+3,j]<-1
      adj.mat[j+3,j+4]<-1
      adj.mat[j+4,j+1]<-1
      adj.mat[j+4,j+2]<-1
      adj.mat[j+4,j+3]<-1
      adj.mat[j+4,j]<-1
      hh.size.vect[j:(j+4)]<-5
    }
    
    temp61<-temp52+1
    temp62<-hh.size.synth$Freq[1]+2*hh.size.synth$Freq[2]+3*hh.size.synth$Freq[3]+4*hh.size.synth$Freq[4]+5*hh.size.synth$Freq[5]+6*hh.size.synth$Freq[6]
    for(j in seq(temp61,temp62,by=6) ){
      adj.mat[j,j+1]<-1
      adj.mat[j,j+2]<-1
      adj.mat[j,j+3]<-1
      adj.mat[j,j+4]<-1
      adj.mat[j,j+5]<-1
      adj.mat[j+1,j]<-1
      adj.mat[j+1,j+2]<-1
      adj.mat[j+1,j+3]<-1
      adj.mat[j+1,j+4]<-1
      adj.mat[j+1,j+5]<-1
      adj.mat[j+2,j+1]<-1
      adj.mat[j+2,j]<-1
      adj.mat[j+2,j+3]<-1
      adj.mat[j+2,j+4]<-1
      adj.mat[j+2,j+5]<-1
      adj.mat[j+3,j+1]<-1
      adj.mat[j+3,j+2]<-1
      adj.mat[j+3,j]<-1
      adj.mat[j+3,j+4]<-1
      adj.mat[j+3,j+5]<-1
      adj.mat[j+4,j+1]<-1
      adj.mat[j+4,j+2]<-1
      adj.mat[j+4,j+3]<-1
      adj.mat[j+4,j]<-1
      adj.mat[j+4,j+5]<-1
      adj.mat[j+5,j+1]<-1
      adj.mat[j+5,j+2]<-1
      adj.mat[j+5,j+3]<-1
      adj.mat[j+5,j+4]<-1
      adj.mat[j+5,j]<-1
      hh.size.vect[j:(j+5)]<-6
    }
    
    temp71<-temp62+1
    temp72<-hh.size.synth$Freq[1]+2*hh.size.synth$Freq[2]+3*hh.size.synth$Freq[3]+4*hh.size.synth$Freq[4]+5*hh.size.synth$Freq[5]+6*hh.size.synth$Freq[6]+7*hh.size.synth$Freq[7]
    for(j in seq(temp71,temp72,by=7) ){
      adj.mat[j,j+1]<-1
      adj.mat[j,j+2]<-1
      adj.mat[j,j+3]<-1
      adj.mat[j,j+4]<-1
      adj.mat[j,j+5]<-1
      adj.mat[j,j+6]<-1
      adj.mat[j+1,j]<-1
      adj.mat[j+1,j+2]<-1
      adj.mat[j+1,j+3]<-1
      adj.mat[j+1,j+4]<-1
      adj.mat[j+1,j+5]<-1
      adj.mat[j+1,j+6]<-1
      adj.mat[j+2,j+1]<-1
      adj.mat[j+2,j]<-1
      adj.mat[j+2,j+3]<-1
      adj.mat[j+2,j+4]<-1
      adj.mat[j+2,j+5]<-1
      adj.mat[j+2,j+6]<-1
      adj.mat[j+3,j+1]<-1
      adj.mat[j+3,j+2]<-1
      adj.mat[j+3,j]<-1
      adj.mat[j+3,j+4]<-1
      adj.mat[j+3,j+5]<-1
      adj.mat[j+3,j+6]<-1
      adj.mat[j+4,j+1]<-1
      adj.mat[j+4,j+2]<-1
      adj.mat[j+4,j+3]<-1
      adj.mat[j+4,j]<-1
      adj.mat[j+4,j+5]<-1
      adj.mat[j+4,j+6]<-1
      adj.mat[j+5,j+1]<-1
      adj.mat[j+5,j+2]<-1
      adj.mat[j+5,j+3]<-1
      adj.mat[j+5,j+4]<-1
      adj.mat[j+5,j]<-1
      adj.mat[j+5,j+6]<-1
      adj.mat[j+6,j+1]<-1
      adj.mat[j+6,j+2]<-1
      adj.mat[j+6,j+3]<-1
      adj.mat[j+6,j+4]<-1
      adj.mat[j+6,j+5]<-1
      adj.mat[j+6,j]<-1
      hh.size.vect[j:(j+6)]<-7
    }
    
    #Constructing hh.id
    hh.id.vect<-NULL
    hh.id.vect<-c(1:hh.size.synth$Freq[1])
    temp.seq<-seq(hh.id.vect[length(hh.id.vect)]+1,hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[2])
    hh.id.vect<-c(hh.id.vect,rep(temp.seq,each=2))
    temp.seq<-seq(hh.id.vect[length(hh.id.vect)]+1,hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[3])
    hh.id.vect<-c(hh.id.vect,rep(temp.seq,each=3))
    temp.seq<-seq(hh.id.vect[length(hh.id.vect)]+1,hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[4])
    hh.id.vect<-c(hh.id.vect,rep(temp.seq,each=4))
    temp.seq<-seq(hh.id.vect[length(hh.id.vect)]+1,hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[5])
    hh.id.vect<-c(hh.id.vect,rep(temp.seq,each=5))
    temp.seq<-seq(hh.id.vect[length(hh.id.vect)]+1,hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[6])
    hh.id.vect<-c(hh.id.vect,rep(temp.seq,each=6))
    temp.seq<-seq(hh.id.vect[length(hh.id.vect)]+1,hh.id.vect[length(hh.id.vect)]+hh.size.synth$Freq[7])
    hh.id.vect<-c(hh.id.vect,rep(temp.seq,each=7))
    vert.attr<-list()
    vert.attr.name<-list()
    vert.attr[[1]]<-hh.size.vect
    vert.attr[[2]]<-hh.id.vect
    vert.attr.name[[1]]<-"hh_size"
    vert.attr.name[[2]]<-"hh_id"
    HH.netw[[w]]<-network(adj.mat, directed = FALSE, vertex.attr = vert.attr, vertex.attrnames =vert.attr.name )
  }
  
  return(HH.netw)
}


