#libraries
library(ggplot2)
library(ggpubr)
library(grid)
######################################################
# Main Text
######################################################
# FIG 1 Het Immunity and pathogen characteristics----
# Het immunity
het.imm<-function(t,t.imm.lim,typeIC){
  if (typeIC==1){
    value<-0.5
  }
  if (typeIC==2){
    value<-(t/t.imm.lim)
  }
  if (typeIC==3){
    value<-0.25+(t/(2*t.imm.lim))
  }
  if (typeIC==4){
    value<-2-t/10
  }
  if (typeIC==5){
      value<-1.5
  }
  if (typeIC==6){
      value<-1.75-t/20
  }
  
  
  return(value)
}
time.points<-seq(0,10,0.1)

hi1<-NULL
hi2<-NULL
hi3<-NULL
hi4<-NULL
hi5<-NULL
hi6<-NULL


for (i in time.points){
  hi1<-c(hi1,het.imm(i,10,1))
  hi2<-c(hi2,het.imm(i,10,2))
  hi3<-c(hi3,het.imm(i,10,3))
}

for (i in time.points){
  hi4<-c(hi4,het.imm(i,10,4))
  hi5<-c(hi5,het.imm(i,10,5))
  hi6<-c(hi6,het.imm(i,10,6))
}


hetImm.ds1<-data.frame(TimePoints=c(rep(time.points,3)), Type=c(rep("Constant",length(time.points)),rep("Decreasing2",length(time.points)),rep("Decreasing1",length(time.points))),value=c(hi1,hi2,hi3))
hetImm.ds4<-data.frame(TimePoints=c(rep(time.points,3)), Type=c(rep("Constant",length(time.points)),rep("Decreasing2",length(time.points)),rep("Decreasing1",length(time.points))),value=c(hi5,hi4,hi6))


library(wesanderson)
pal2<- wes_palette("Chevalier1",4,type = "discrete")
pal<-c(pal2[1],"steelblue",pal2[3],pal2[4])


v1<-ggplot(hetImm.ds1, aes(x=TimePoints, y=value, col=Type))+geom_line( lwd=0.9,lty=1)+scale_color_manual(values=pal)+
  theme_classic()+theme(
    legend.key.size = unit(0.8, "cm"),
    legend.key.width = unit(0.8,"cm"),
    legend.title = element_text(size=14, face = "bold"),
    legend.text = element_text(size=13),
    #    legend.position = "top",
    axis.text.x = element_text(size=13),
    axis.title = element_text(size=14),
    axis.text.y = element_text(size=13),
    text = element_text(size=17, face = "bold"),
  )+xlab("Days since recovery")+ylab("Susceptibility ")+ylim(c(0,1))+labs(color="CompEff:")


v2<-ggplot(hetImm.ds4, aes(x=TimePoints, y=value, col=Type))+geom_line( lwd=0.9,lty=1)+scale_color_brewer(palette="Dark2")+
  theme_classic()+theme(
    legend.key.size = unit(0.8, "cm"),
    legend.key.width = unit(0.8,"cm"),
    legend.title = element_text(size=14, face = "bold"),
    legend.text = element_text(size=13),
    #    legend.position = "top",
    axis.text.x = element_text(size=13),
    axis.title = element_text(size=14),
    axis.text.y = element_text(size=13),
    text = element_text(size=17, face = "bold"),
  )+xlab("Days since recovery")+ylab("Susceptibility ")+ylim(c(1,2))+labs(color="CoopEff:")

plot(v1)
ggarrange(v1,v2, labels = c("A","B"))
#saved: 1000*500







# FIG 2 Independence + Home isolation ----
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(24203)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t254_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","13","26","39"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","13","26","39"),z=c(0,18,36,54))

aa2<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Beginning","13"="Ascending","26"="Peak","39"="Descending")))
plot(aa2)

# FLU-COVID - Compliance: COVID 1, FLU 0 - Prevalence curves 
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(123)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","13","26","39"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","13","26","39"),z=c(0,25,50,75))

aa21<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Beginning","13"="Ascending","26"="Peak","39"="Descending")))+ylim(c(0,0.25))

ggarrange(aa2,aa21,
          nrow=2,
          labels = c("A","B"), common.legend = TRUE)


#Saved: 950x820


# FIG 3 Home Iso + Het immunity (attack rate + cases at peak) ----
# FLU-COVID - Compliance: COVID 1, FLU 0 - Baseline 2 
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.Fig3a.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.Fig3a.AR$Pathogen<-factor(dt.Fig3a.AR$Pathogen, levels = c("Path1","Path2"))
dt.Fig3a.AR$IntroductionTime<-factor(dt.Fig3a.AR$IntroductionTime, levels = c("0","13","26","39"))


aa101<-ggplot(data = dt.Fig3a.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.7))

temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(123)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Fig3a.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.Fig3a.CP$Pathogen<-factor(dt.Fig3a.CP$Pathogen, levels = c("Path1","Path2"))
dt.Fig3a.CP$IntroductionTime<-factor(dt.Fig3a.CP$IntroductionTime, levels = c("0","13","26","39"))


aa31<-ggplot(data = dt.Fig3a.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","FLU"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ylim(c(0,0.25))

# FLU-COVID - Cross-immunity: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.bc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/18/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/27/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.Fig3b.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.Fig3b.AR$Pathogen<-factor(dt.Fig3b.AR$Pathogen, levels = c("Path1","Path2"))
dt.Fig3b.AR$IntroductionTime<-factor(dt.Fig3b.AR$IntroductionTime, levels = c("0","13","26","39"))


aa10352<-ggplot(data = dt.Fig3b.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.7))


# FLU-COVID - Prevalence curves over time
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(123)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/18/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/27/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)


dt.Fig3b.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.Fig3b.CP$Pathogen<-factor(dt.Fig3b.CP$Pathogen, levels = c("Path1","Path2"))
dt.Fig3b.CP$IntroductionTime<-factor(dt.Fig3b.CP$IntroductionTime, levels = c("0","13","26","39"))


aa32<-ggplot(data = dt.Fig3b.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ylim(c(0,0.25))

ggarrange(aa101,aa10352,aa31,aa32, ncol = 2,nrow = 2, labels = c("A", "B", "C","D"),common.legend = TRUE)


#Table
#Figure 1
print("Panel-A")
tapply(dt.Fig3a.AR$FinalSize/2500, INDEX = list(dt.Fig3a.AR$Pathogen,dt.Fig3a.AR$IntroductionTime),FUN = mean)
tapply(dt.Fig3a.AR$FinalSize/2500, INDEX = list(dt.Fig3a.AR$Pathogen,dt.Fig3a.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig3a.AR$FinalSize/2500, INDEX = list(dt.Fig3a.AR$Pathogen,dt.Fig3a.AR$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 2
print("Panel-B")
tapply(dt.Fig3b.AR$FinalSize/2500, INDEX = list(dt.Fig3b.AR$Pathogen,dt.Fig3b.AR$IntroductionTime),FUN = mean)
tapply(dt.Fig3b.AR$FinalSize/2500, INDEX = list(dt.Fig3b.AR$Pathogen,dt.Fig3b.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig3b.AR$FinalSize/2500, INDEX = list(dt.Fig3b.AR$Pathogen,dt.Fig3b.AR$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 3
print("Panel-C")
tapply(dt.Fig3a.CP$PeakCases/2500, INDEX = list(dt.Fig3a.CP$Pathogen,dt.Fig3a.CP$IntroductionTime),FUN = mean)
tapply(dt.Fig3a.CP$PeakCases/2500, INDEX = list(dt.Fig3a.CP$Pathogen,dt.Fig3a.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig3a.CP$PeakCases/2500, INDEX = list(dt.Fig3a.CP$Pathogen,dt.Fig3a.CP$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 4
print("Panel-D")
tapply(dt.Fig3b.CP$PeakCases/2500, INDEX = list(dt.Fig3b.CP$Pathogen,dt.Fig3b.CP$IntroductionTime),FUN = mean)
tapply(dt.Fig3b.CP$PeakCases/2500, INDEX = list(dt.Fig3b.CP$Pathogen,dt.Fig3b.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig3b.CP$PeakCases/2500, INDEX = list(dt.Fig3b.CP$Pathogen,dt.Fig3b.CP$IntroductionTime),FUN = quantile, p=0.025 )

#Saved: 950x820

# FIG 4 Long-term Competition + Cooperation ----
temp.path<-NULL
temp.bc<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/19/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/28/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.Fig4a.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.Fig4a.AR$Pathogen<-factor(dt.Fig4a.AR$Pathogen, levels = c("Path1","Path2"))
dt.Fig4a.AR$IntroductionTime<-factor(dt.Fig4a.AR$IntroductionTime, levels = c("0","13","26","39"))


aa101<-ggplot(data = dt.Fig4a.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.9))


# FLU-COVID - Prevalence curves over time
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(12345)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/19/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/28/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","13","26","39"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","13","26","39"),z=c(0,25,50,75))

aa21<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","FLU"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "grey",
                                  size = 0.2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=10),
  legend.position = "top",
  legend.title = element_text(size=11, face = "bold"),
  axis.text.x = element_text(size=10, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=11),
  axis.text.y = element_text(size=10),
  strip.text = element_text(size = 11)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Beginning","13"="Ascending","26"="Peak","39"="Descending")))+ylim(c(0,0.3))

dt.Fig4a.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.Fig4a.CP$Pathogen<-factor(dt.Fig4a.CP$Pathogen, levels = c("Path1","Path2"))
dt.Fig4a.CP$IntroductionTime<-factor(dt.Fig4a.CP$IntroductionTime, levels = c("0","13","26","39"))


aa31<-ggplot(data = dt.Fig4a.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ ylim(c(0,0.25))

# FLU-COVID -  2 co-op 

temp.path<-NULL
temp.bc<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/17/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/26/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/35/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.Fig4b.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.bc)
dt.Fig4b.AR$Pathogen<-factor(dt.Fig4b.AR$Pathogen, levels = c("Path2","Path1"))
dt.Fig4b.AR$IntroductionTime<-factor(dt.Fig4b.AR$IntroductionTime, levels = c("0","13","26","39"))


aa102<-ggplot(data = dt.Fig4b.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#E69F00","#56B4E9"), labels=c("Influenza","COVID-19"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.9))


# FLU-COVID - Prevalence curves over time
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(123)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/17/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/26/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompSM/35/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path2","Path1"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","13","26","39"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","13","26","39"),z=c(0,12.5,25,37.5))

aa22<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","FLU"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "grey",
                                  size = 0.2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Beginning","13"="Ascending","26"="Peak","39"="Descending")))+ylim(c(0,0.3))

dt.Fig4b.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.Fig4b.CP$Pathogen<-factor(dt.Fig4b.CP$Pathogen, levels = c("Path2","Path1"))
dt.Fig4b.CP$IntroductionTime<-factor(dt.Fig4b.CP$IntroductionTime, levels = c("0","13","26","39"))


aa32<-ggplot(data = dt.Fig4b.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#E69F00","#56B4E9"), labels=c("Influenza","COVID-19"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ylim(c(0,0.25))

ggarrange(aa101,aa102,aa31,aa32, nrow = 2, ncol = 2, labels = c("A","B","C","D"), common.legend = TRUE)

#Saved: 950 x 820
#Table
#Figure 1
print("A")
tapply(dt.Fig4a.AR$FinalSize/2500, INDEX = list(dt.Fig4a.AR$Pathogen,dt.Fig4a.AR$IntroductionTime),FUN = mean)
tapply(dt.Fig4a.AR$FinalSize/2500, INDEX = list(dt.Fig4a.AR$Pathogen,dt.Fig4a.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig4a.AR$FinalSize/2500, INDEX = list(dt.Fig4a.AR$Pathogen,dt.Fig4a.AR$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 2
print("B")
tapply(dt.Fig4b.AR$FinalSize/2500, INDEX = list(dt.Fig4b.AR$Pathogen,dt.Fig4b.AR$IntroductionTime),FUN = mean)
tapply(dt.Fig4b.AR$FinalSize/2500, INDEX = list(dt.Fig4b.AR$Pathogen,dt.Fig4b.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig4b.AR$FinalSize/2500, INDEX = list(dt.Fig4b.AR$Pathogen,dt.Fig4b.AR$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 4
print("C")
tapply(dt.Fig4a.CP$PeakCases/2500, INDEX = list(dt.Fig4a.CP$Pathogen,dt.Fig4a.CP$IntroductionTime),FUN = mean)
tapply(dt.Fig4a.CP$PeakCases/2500, INDEX = list(dt.Fig4a.CP$Pathogen,dt.Fig4a.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig4a.CP$PeakCases/2500, INDEX = list(dt.Fig4a.CP$Pathogen,dt.Fig4a.CP$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 4
print("D")
tapply(dt.Fig4b.CP$PeakCases/2500, INDEX = list(dt.Fig4b.CP$Pathogen,dt.Fig4b.CP$IntroductionTime),FUN = mean)
tapply(dt.Fig4b.CP$PeakCases/2500, INDEX = list(dt.Fig4b.CP$Pathogen,dt.Fig4b.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.Fig4b.CP$PeakCases/2500, INDEX = list(dt.Fig4b.CP$Pathogen,dt.Fig4b.CP$IntroductionTime),FUN = quantile, p=0.025 )


######################################################
# Independence + Home isolation ----
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(24203)
temp.pk<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path2",length(epi.det$Days)),rep("Path1",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("1",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("2",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("3",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("4",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)



dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","1","2","3","4"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","1","2","3","4"),z=c(0,0,18,36,108))

IND<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
#  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
#  strip.text.x = element_text(size=15, face = "bold"),
  panel.spacing = unit(0.75, "cm", data = NULL),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
#  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Precedent","1"="Beginning","2"="Ascending","3"="Peak","4"="Subsequent")))
plot(IND)

# FLU-COVID - COVID 1, FLU 0 - Prevalence curves 
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(24203)
temp.pk<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path2",length(epi.det$Days)),rep("Path1",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("1",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("2",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("3",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("4",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)



dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","1","2","3","4"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","1","2","3","4"),z=c(0,0,18,36,108))

SCI<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Precedent","1"="Beginning","2"="Ascending","3"="Peak","4"="Subsequent")))
#plot(SCI)

ggarrange(IND,SCI,
          nrow=2,
          labels = c("A","B"), common.legend = TRUE)


#Saved: 950x820
# Independence, STI and LTI

temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.sce<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("3",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("3",2*length(epi.outbreak)))





load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("3",2*length(epi.outbreak)))





load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("3",2*length(epi.outbreak)))




load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce,rep("3",2*length(epi.outbreak)))


dt.BSLINT.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Scenario"=temp.sce)
dt.BSLINT.AR$Pathogen<-factor(dt.BSLINT.AR$Pathogen, levels = c("Path1","Path2"))
dt.BSLINT.AR$IntroductionTime<-factor(dt.BSLINT.AR$IntroductionTime, levels = c("0","1","2","3","4"))


BSLINT.AR<-ggplot(data = dt.BSLINT.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
  strip.text.x = element_text(size=15, face = "bold"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15),
  panel.spacing = unit(0.75, "cm", data = NULL),
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.91))+facet_grid(~Scenario,  scales = "fixed", labeller = labeller(Scenario=c("0"="Independent","1"="Short-term cross-protection","2"="Long-term competition","3"="Long-term cooperation")))

round(tapply(dt.BSLINT.AR$FinalSize/2500, INDEX = list(dt.BSLINT.AR$Pathogen,dt.BSLINT.AR$IntroductionTime, dt.BSLINT.AR$Scenario),FUN = mean),2)
round(tapply(dt.BSLINT.AR$FinalSize/2500, INDEX = list(dt.BSLINT.AR$Pathogen,dt.BSLINT.AR$IntroductionTime, dt.BSLINT.AR$Scenario),FUN = quantile, p=0.975 ),2)
round(tapply(dt.BSLINT.AR$FinalSize/2500, INDEX = list(dt.BSLINT.AR$Pathogen,dt.BSLINT.AR$IntroductionTime, dt.BSLINT.AR$Scenario),FUN = quantile, p=0.025 ),2)


plot(BSLINT.AR)

temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL
temp.sce.pk<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk,rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk,rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk,rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk,rep("3",2))
}


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk,rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk,rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk,rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk,rep("3",2))
}





load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk,rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk,rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk,rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk,rep("3",2))
}





load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk,rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk,rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk,rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk,rep("3",2))
}




load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk,rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk,rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOMP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk,rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTLCOOP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk,rep("3",2))
}


dt.BSLINT.CP<-data.frame("FinalSize"=temp.pk,"Pathogen"=temp.path.pk,"IntroductionTime"=temp.t2.pk, "Scenario"=temp.sce.pk)
dt.BSLINT.CP$Pathogen<-factor(dt.BSLINT.CP$Pathogen, levels = c("Path1","Path2"))
dt.BSLINT.CP$IntroductionTime<-factor(dt.BSLINT.CP$IntroductionTime, levels = c("0","1","2","3","4"))


BSLINT.CP<-ggplot(data = dt.BSLINT.CP, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
  strip.text.x = element_text(size=15, face = "bold"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15),
  panel.spacing = unit(0.75, "cm", data = NULL),
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.45))+facet_grid(~Scenario,  scales = "fixed", labeller = labeller(Scenario=c("0"="Independent","1"="Short-term cross-protection","2"="Long-term competition","3"="Long-term cooperation")))

plot(BSLINT.CP)

round(tapply(dt.BSLINT.CP$FinalSize/2500, INDEX = list(dt.BSLINT.CP$Pathogen,dt.BSLINT.CP$IntroductionTime, dt.BSLINT.CP$Scenario),FUN = mean),2)
round(tapply(dt.BSLINT.CP$FinalSize/2500, INDEX = list(dt.BSLINT.CP$Pathogen,dt.BSLINT.CP$IntroductionTime, dt.BSLINT.CP$Scenario),FUN = quantile, p=0.975 ),2)
round(tapply(dt.BSLINT.CP$FinalSize/2500, INDEX = list(dt.BSLINT.CP$Pathogen,dt.BSLINT.CP$IntroductionTime, dt.BSLINT.CP$Scenario),FUN = quantile, p=0.025 ),2)






# # Home Iso: Independent, short-term competition and long term competition or cooperation (attack rate + cases at peak) ----
# FLU-COVID - Independent  
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.sce<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("0",2*length(epi.outbreak)))


dt.FigNHIIND.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHIIND.AR$Pathogen<-factor(dt.FigNHIIND.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHIIND.AR$IntroductionTime<-factor(dt.FigNHIIND.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.IND.AR<-ggplot(data = dt.FigNHIIND.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))

plot(HI.IND.AR)
round(tapply(dt.FigNHIIND.AR$FinalSize/2500, INDEX = list(dt.FigNHIIND.AR$Pathogen,dt.FigNHIIND.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigNHIIND.AR$FinalSize/2500, INDEX = list(dt.FigNHIIND.AR$Pathogen,dt.FigNHIIND.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)
round(tapply(dt.FigNHIIND.AR$FinalSize/2500, INDEX = list(dt.FigNHIIND.AR$Pathogen,dt.FigNHIIND.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL
temp.sce.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk, rep("0",2))
}


dt.FigNHIIND.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHIIND.CP$Pathogen<-factor(dt.FigNHIIND.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHIIND.CP$IntroductionTime<-factor(dt.FigNHIIND.CP$IntroductionTime, levels = c("0","1","2","3","4"))


HI.IND.CP<-ggplot(data = dt.FigNHIIND.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.25))
plot(HI.IND.CP)

# FLU-COVID - Short-term cross-immunity: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))


dt.FigNHSTC.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHSTC.AR$Pathogen<-factor(dt.FigNHSTC.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHSTC.AR$IntroductionTime<-factor(dt.FigNHSTC.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.STC.AR<-ggplot(data = dt.FigNHSTC.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)
plot(HI.STC.AR)

temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}


dt.FigNHISTC.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHISTC.CP$Pathogen<-factor(dt.FigNHISTC.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHISTC.CP$IntroductionTime<-factor(dt.FigNHISTC.CP$IntroductionTime, levels = c("0","1","2","3","4"))

HI.STC.CP<-ggplot(data = dt.FigNHISTC.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.25))
plot(HI.STC.CP)

# FLU-COVID - Long-term competition: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("2",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("2",2*length(epi.outbreak)))


dt.FigNHLTComp.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHLTComp.AR$Pathogen<-factor(dt.FigNHLTComp.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHLTComp.AR$IntroductionTime<-factor(dt.FigNHLTComp.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.LTComp.AR<-ggplot(data = dt.FigNHLTComp.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)
plot(HI.LTComp.AR)

tapply(dt.FigNHLTComp.AR$FinalSize/2500, INDEX = list(dt.FigNHLTComp.AR$Pathogen,dt.FigNHLTComp.AR$IntroductionTime),FUN = mean)
tapply(dt.FigNHLTComp.AR$FinalSize/2500, INDEX = list(dt.FigNHLTComp.AR$Pathogen,dt.FigNHLTComp.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigNHLTComp.AR$FinalSize/2500, INDEX = list(dt.FigNHLTComp.AR$Pathogen,dt.FigNHLTComp.AR$IntroductionTime),FUN = quantile, p=0.025 )


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOMP/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk, rep("2",2))
}


dt.FigNHILTComp.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHILTComp.CP$Pathogen<-factor(dt.FigNHILTComp.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHILTComp.CP$IntroductionTime<-factor(dt.FigNHILTComp.CP$IntroductionTime, levels = c("0","1","2","3","4"))


HI.LTComp.CP<-ggplot(data = dt.FigNHILTComp.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.25))
plot(HI.LTComp.CP)


tapply(dt.FigNHILTComp.CP$PeakCases/2500, INDEX = list(dt.FigNHILTComp.CP$Pathogen,dt.FigNHILTComp.CP$IntroductionTime),FUN = mean)
tapply(dt.FigNHILTComp.CP$PeakCases/2500, INDEX = list(dt.FigNHILTComp.CP$Pathogen,dt.FigNHILTComp.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigNHILTComp.CP$PeakCases/2500, INDEX = list(dt.FigNHILTComp.CP$Pathogen,dt.FigNHILTComp.CP$IntroductionTime),FUN = quantile, p=0.025 )

# FLU-COVID - Long-term co-operation: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("3",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("3",2*length(epi.outbreak)))


dt.FigNHLTCoop.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHLTCoop.AR$Pathogen<-factor(dt.FigNHLTCoop.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHLTCoop.AR$IntroductionTime<-factor(dt.FigNHLTCoop.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.LTCoop.AR<-ggplot(data = dt.FigNHLTCoop.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)
plot(HI.LTCoop.AR)

tapply(dt.FigNHLTCoop.AR$FinalSize/2500, INDEX = list(dt.FigNHLTCoop.AR$Pathogen,dt.FigNHLTCoop.AR$IntroductionTime),FUN = mean)
tapply(dt.FigNHLTCoop.AR$FinalSize/2500, INDEX = list(dt.FigNHLTCoop.AR$Pathogen,dt.FigNHLTCoop.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigNHLTCoop.AR$FinalSize/2500, INDEX = list(dt.FigNHLTCoop.AR$Pathogen,dt.FigNHLTCoop.AR$IntroductionTime),FUN = quantile, p=0.025 )


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOCOOP/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk, rep("3",2))
}


dt.FigNHILTCoop.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHILTCoop.CP$Pathogen<-factor(dt.FigNHILTCoop.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHILTCoop.CP$IntroductionTime<-factor(dt.FigNHILTCoop.CP$IntroductionTime, levels = c("0","1","2","3","4"))


HI.LTCoop.CP<-ggplot(data = dt.FigNHILTCoop.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.25))

plot(HI.LTCoop.CP)

tapply(dt.FigNHILTCoop.CP$PeakCases/2500, INDEX = list(dt.FigNHILTCoop.CP$Pathogen,dt.FigNHILTCoop.CP$IntroductionTime),FUN = mean)
tapply(dt.FigNHILTCoop.CP$PeakCases/2500, INDEX = list(dt.FigNHILTCoop.CP$Pathogen,dt.FigNHILTCoop.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigNHILTCoop.CP$PeakCases/2500, INDEX = list(dt.FigNHILTCoop.CP$Pathogen,dt.FigNHILTCoop.CP$IntroductionTime),FUN = quantile, p=0.025 )

#Saved: 950x820


# Plotting together


dt.HIALL.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Scenario"=temp.sce)
dt.HIALL.AR$Pathogen<-factor(dt.HIALL.AR$Pathogen, levels = c("Path1","Path2"))
dt.HIALL.AR$IntroductionTime<-factor(dt.HIALL.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HIALL.AR<-ggplot(data = dt.HIALL.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
  strip.text.x = element_text(size=15, face = "bold"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15),
  panel.spacing = unit(0.75, "cm", data = NULL),
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))+facet_grid(~Scenario,  scales = "fixed", labeller = labeller(Scenario=c("0"="Independent","1"="Short-term cross-protection","2"="Long-term competition","3"="Long-term cooperation")))
plot(HIALL.AR)

round(tapply(dt.HIALL.AR$FinalSize/2500, INDEX = list(dt.HIALL.AR$Pathogen,dt.HIALL.AR$IntroductionTime, dt.HIALL.AR$Scenario),FUN = mean),2)
round(tapply(dt.HIALL.AR$FinalSize/2500, INDEX = list(dt.HIALL.AR$Pathogen,dt.HIALL.AR$IntroductionTime, dt.HIALL.AR$Scenario),FUN = quantile, p=0.975 ),2)
round(tapply(dt.HIALL.AR$FinalSize/2500, INDEX = list(dt.HIALL.AR$Pathogen,dt.HIALL.AR$IntroductionTime, dt.HIALL.AR$Scenario),FUN = quantile, p=0.025 ),2)



dt.HIALL.CP<-data.frame("FinalSize"=temp.pk,"Pathogen"=temp.path.pk,"IntroductionTime"=temp.t2.pk, "Scenario"=temp.sce.pk)
dt.HIALL.CP$Pathogen<-factor(dt.HIALL.CP$Pathogen, levels = c("Path1","Path2"))
dt.HIALL.CP$IntroductionTime<-factor(dt.HIALL.CP$IntroductionTime, levels = c("0","1","2","3","4"))


HIALL.CP<-ggplot(data = dt.HIALL.CP, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
  strip.text.x = element_text(size=15, face = "bold"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15),
  panel.spacing = unit(0.75, "cm", data = NULL),
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.23))+facet_grid(~Scenario,  scales = "fixed", labeller = labeller(Scenario=c("0"="Independent","1"="Short-term cross-protection","2"="Long-term competition","3"="Long-term cooperation")))
plot(HIALL.CP)


# FLU-COVID - Independent - BC influenza  
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigNHIIND.BC.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHIIND.BC.AR$Pathogen<-factor(dt.FigNHIIND.BC.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHIIND.BC.AR$IntroductionTime<-factor(dt.FigNHIIND.BC.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.IND.BC.AR<-ggplot(data = dt.FigNHIIND.BC.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))+geom_hline(yintercept = 0.1, lty=2, alpha=0.5)
plot(HI.IND.BC.AR)


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigNHIIND.BC.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHIIND.BC.CP$Pathogen<-factor(dt.FigNHIIND.BC.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHIIND.BC.CP$IntroductionTime<-factor(dt.FigNHIIND.BC.CP$IntroductionTime, levels = c("0","1","2","3","4"))


HI.IND.BC.CP<-ggplot(data = dt.FigNHIIND.BC.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.25))
plot(HI.IND.BC.CP)

temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(24203)
temp.pk<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path2",length(epi.det$Days)),rep("Path1",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("1",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("2",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("3",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOIND/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("4",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)



dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","1","2","3","4"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","1","2","3","4"),z=c(0,0,25,50,125))

HIFLU.Prev<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Days")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Precedent","1"="Beginning","2"="Ascending","3"="Peak","4"="Subsequent")))

plot(HIFLU.Prev)
ggarrange(HIFLU.Prev,ggarrange(HI.IND.BC.AR,HI.IND.BC.CP, labels = c("B,C"), common.legend = TRUE), labels = c("A","B","C"), common.legend = TRUE, nrow = 2, heights = c(0.8,1.2))
#









# # Global contact reduction (TODO): Independent, short-term competition (attack rate) ----
# FLU-COVID - Independent  
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.cd<-NULL
temp.f1<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))

temp.f1<-c(temp.f1,rep("0",length(temp.fs)))



dt.FigGCRIND.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Reduction"=temp.cd)
dt.FigGCRIND.AR$Pathogen<-factor(dt.FigGCRIND.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigGCRIND.AR$IntroductionTime<-factor(dt.FigGCRIND.AR$IntroductionTime, levels = c("0","1","2","3","4"))
dt.FigGCRIND.AR$Reduction<-factor(dt.FigGCRIND.AR$Reduction, levels = c("2","0","1"))


GCR.IND.AR<-ggplot(data = dt.FigGCRIND.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.91))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30% reduction","1"="65% reduction")))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)
plot(GCR.IND.AR)



round(tapply(dt.FigGCRIND.AR$FinalSize/2500, INDEX = list(dt.FigGCRIND.AR$Pathogen,dt.FigGCRIND.AR$IntroductionTime, dt.FigGCRIND.AR$Reduction),FUN = mean),2)
round(tapply(dt.FigGCRIND.AR$FinalSize/2500, INDEX = list(dt.FigGCRIND.AR$Pathogen,dt.FigGCRIND.AR$IntroductionTime, dt.FigGCRIND.AR$Reduction),FUN = quantile, p=0.975 ),2)
round(tapply(dt.FigGCRIND.AR$FinalSize/2500, INDEX = list(dt.FigGCRIND.AR$Pathogen,dt.FigGCRIND.AR$IntroductionTime, dt.FigGCRIND.AR$Reduction),FUN = quantile, p=0.025 ),2)


# FLU-COVID - Cross-immunity during infection  
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.cd<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/NOINTSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2108_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



dt.FigGCRSCI.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Reduction"=temp.cd)
dt.FigGCRSCI.AR$Pathogen<-factor(dt.FigGCRSCI.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigGCRSCI.AR$IntroductionTime<-factor(dt.FigGCRSCI.AR$IntroductionTime, levels = c("0","1","2","3","4"))
dt.FigGCRSCI.AR$Reduction<-factor(dt.FigGCRSCI.AR$Reduction, levels = c("2","0","1"))


GCR.SCI.AR<-ggplot(data = dt.FigGCRSCI.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.91))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30% reduction","1"="65% reduction")))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)
plot(GCR.SCI.AR)

round(tapply(dt.FigGCRSCI.AR$FinalSize/2500, INDEX = list(dt.FigGCRSCI.AR$Pathogen,dt.FigGCRSCI.AR$IntroductionTime, dt.FigGCRSCI.AR$Reduction),FUN = mean),2)
round(tapply(dt.FigGCRSCI.AR$FinalSize/2500, INDEX = list(dt.FigGCRSCI.AR$Pathogen,dt.FigGCRSCI.AR$IntroductionTime, dt.FigGCRSCI.AR$Reduction),FUN = quantile, p=0.975 ),2)
round(tapply(dt.FigGCRSCI.AR$FinalSize/2500, INDEX = list(dt.FigGCRSCI.AR$Pathogen,dt.FigGCRSCI.AR$IntroductionTime, dt.FigGCRSCI.AR$Reduction),FUN = quantile, p=0.025 ),2)


#
# 
# # 


require(grid)   # for the textGrob() function

GCR.AR <- ggarrange(GCR.IND.AR + rremove("xlab") + rremove("ylab"), GCR.SCI.AR + rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                    labels = c("A","B"),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "top",
                    align = "hv", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(GCR.AR, left = textGrob("Attack Rate", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))





# # Same R*: Independent, short-term competition and long term competition or cooperation (attack rate + cases at peak) ----
# FLU-COVID - Independent
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARIND.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARIND.AR$Pathogen<-factor(dt.FigRSTARIND.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARIND.AR$IntroductionTime<-factor(dt.FigRSTARIND.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.IND.AR<-ggplot(data = dt.FigRSTARIND.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
  strip.text.x = element_text(size=15, face = "bold"),
  panel.spacing = unit(0.75, "cm", data = NULL),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.6))
plot(RSTAR.IND.AR)

round(tapply(dt.FigRSTARIND.AR$FinalSize/2500, INDEX = list(dt.FigRSTARIND.AR$Pathogen,dt.FigRSTARIND.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARIND.AR$FinalSize/2500, INDEX = list(dt.FigRSTARIND.AR$Pathogen,dt.FigRSTARIND.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARIND.AR$FinalSize/2500, INDEX = list(dt.FigRSTARIND.AR$Pathogen,dt.FigRSTARIND.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)



temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARIND/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARIND.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARIND.CP$Pathogen<-factor(dt.FigRSTARIND.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARIND.CP$IntroductionTime<-factor(dt.FigRSTARIND.CP$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.IND.CP<-ggplot(data = dt.FigRSTARIND.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  strip.background = element_rect(color="white", fill="white", size=0.1, linetype="solid"),
  strip.text.x = element_text(size=15, face = "bold"),
  panel.spacing = unit(0.75, "cm", data = NULL),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.2))
plot(RSTAR.IND.CP)

RSTAR.IND<- ggarrange(RSTAR.IND.AR + rremove("xlab"), RSTAR.IND.CP + rremove("xlab"), # remove axis labels from plots
                    labels = c("A","B"),
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "top",
                    align = "hv", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(RSTAR.IND,
                bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))




# FLU-COVID - Short-term cross-immunity: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARSTC.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARSTC.AR$Pathogen<-factor(dt.FigRSTARSTC.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARSTC.AR$IntroductionTime<-factor(dt.FigRSTARSTC.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.STC.AR<-ggplot(data = dt.FigRSTARSTC.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.6))
plot(RSTAR.STC.AR)

round(tapply(dt.FigRSTARSTC.AR$FinalSize/2500, INDEX = list(dt.FigRSTARSTC.AR$Pathogen,dt.FigRSTARSTC.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARSTC.AR$FinalSize/2500, INDEX = list(dt.FigRSTARSTC.AR$Pathogen,dt.FigRSTARSTC.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARSTC.AR$FinalSize/2500, INDEX = list(dt.FigRSTARSTC.AR$Pathogen,dt.FigRSTARSTC.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)



temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARSTC.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARSTC.CP$Pathogen<-factor(dt.FigRSTARSTC.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARSTC.CP$IntroductionTime<-factor(dt.FigRSTARSTC.CP$IntroductionTime, levels = c("0","1","2","3","4"))

RSTAR.STC.CP<-ggplot(data = dt.FigRSTARSTC.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.2))
plot(RSTAR.STC.CP)

# FLU-COVID - Long-term competition: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARLTComp.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARLTComp.AR$Pathogen<-factor(dt.FigRSTARLTComp.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp.AR$IntroductionTime<-factor(dt.FigRSTARLTComp.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp.AR<-ggplot(data = dt.FigRSTARLTComp.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.6))
plot(RSTAR.LTComp.AR)

round(tapply(dt.FigRSTARLTComp.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.AR$Pathogen,dt.FigRSTARLTComp.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARLTComp.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.AR$Pathogen,dt.FigRSTARLTComp.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARLTComp.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.AR$Pathogen,dt.FigRSTARLTComp.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)



temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARLTComp.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARLTComp.CP$Pathogen<-factor(dt.FigRSTARLTComp.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp.CP$IntroductionTime<-factor(dt.FigRSTARLTComp.CP$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp.CP<-ggplot(data = dt.FigRSTARLTComp.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.2))

plot(RSTAR.LTComp.CP)

tapply(dt.FigRSTARLTComp.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp.CP$Pathogen,dt.FigRSTARLTComp.CP$IntroductionTime),FUN = mean)
tapply(dt.FigRSTARLTComp.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp.CP$Pathogen,dt.FigRSTARLTComp.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigRSTARLTComp.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp.CP$Pathogen,dt.FigRSTARLTComp.CP$IntroductionTime),FUN = quantile, p=0.025 )



# FLU-COVID - Long-term competition: COVID 1, FLU 0 + HI - LLI 30-days
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARLTComp30.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARLTComp30.AR$Pathogen<-factor(dt.FigRSTARLTComp30.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp30.AR$IntroductionTime<-factor(dt.FigRSTARLTComp30.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp30.AR<-ggplot(data = dt.FigRSTARLTComp30.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.6))
plot(RSTAR.LTComp30.AR)


round(tapply(dt.FigRSTARLTComp30.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp30.AR$Pathogen,dt.FigRSTARLTComp30.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARLTComp30.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp30.AR$Pathogen,dt.FigRSTARLTComp30.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARLTComp30.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp30.AR$Pathogen,dt.FigRSTARLTComp30.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARLTComp30.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARLTComp30.CP$Pathogen<-factor(dt.FigRSTARLTComp30.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp30.CP$IntroductionTime<-factor(dt.FigRSTARLTComp30.CP$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp30.CP<-ggplot(data = dt.FigRSTARLTComp30.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.2))

tapply(dt.FigRSTARLTComp30.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp30.CP$Pathogen,dt.FigRSTARLTComp30.CP$IntroductionTime),FUN = mean)
tapply(dt.FigRSTARLTComp30.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp30.CP$Pathogen,dt.FigRSTARLTComp30.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigRSTARLTComp30.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp30.CP$Pathogen,dt.FigRSTARLTComp30.CP$IntroductionTime),FUN = quantile, p=0.025 )

RSTAR.COMP.AR <- ggarrange(RSTAR.IND.AR + rremove("xlab")+ rremove("ylab"), RSTAR.STC.AR + rremove("ylab") + rremove("ylab")+ rremove("xlab"),RSTAR.LTComp.AR + rremove("xlab") + rremove("ylab"),RSTAR.LTComp30.AR + rremove("xlab") + rremove("ylab"), # remove axis labels from plots
                    labels = c("Independent","Shor-term cross protection","Long-term competition - baseline","Long-term competition - high competition"),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "top",
                    align = "hv", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(RSTAR.COMP.AR, left = textGrob("Attack Rate", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))

RSTAR.COMP.CP <- ggarrange(RSTAR.IND.CP + rremove("xlab") + rremove("ylab"), RSTAR.STC.CP + rremove("ylab") + rremove("xlab"),RSTAR.LTComp.CP + rremove("xlab") + rremove("ylab"),RSTAR.LTComp30.CP + rremove("xlab") + rremove("ylab"), # remove axis labels from plots
                           labels = c("Independent","Shor-term cross protection","Long-term competition - baseline","Long-term competition - high competition"),
                           ncol = 2, nrow = 2,
                           common.legend = TRUE, legend = "none",
                           align = "hv", 
                           font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(RSTAR.COMP.CP, left = textGrob("Cases at Peak", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))

ggarrange(RSTAR.COMP.AR,RSTAR.COMP.CP, common.legend = TRUE, ncol = 1, legend = "top")





# FLU-COVID - Long-term competition: COVID 1, FLU 0 + HI - R3.3
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP233/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARLTComp.R33.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARLTComp.R33.AR$Pathogen<-factor(dt.FigRSTARLTComp.R33.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp.R33.AR$IntroductionTime<-factor(dt.FigRSTARLTComp.R33.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp.R33.AR<-ggplot(data = dt.FigRSTARLTComp.R33.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0.65,0.92))
plot(RSTAR.LTComp.R33.AR)


round(tapply(dt.FigRSTARLTComp.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.R33.AR$Pathogen,dt.FigRSTARLTComp.R33.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARLTComp.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.R33.AR$Pathogen,dt.FigRSTARLTComp.R33.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARLTComp.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.R33.AR$Pathogen,dt.FigRSTARLTComp.R33.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP233/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARLTComp.R33.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARLTComp.R33.CP$Pathogen<-factor(dt.FigRSTARLTComp.R33.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp.R33.CP$IntroductionTime<-factor(dt.FigRSTARLTComp.R33.CP$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp.R33.CP<-ggplot(data = dt.FigRSTARLTComp.R33.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0.08,0.65))

plot(RSTAR.LTComp.R33.CP)

round(tapply(dt.FigRSTARLTComp.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.R33.AR$Pathogen,dt.FigRSTARLTComp.R33.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARLTComp.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.R33.AR$Pathogen,dt.FigRSTARLTComp.R33.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARLTComp.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp.R33.AR$Pathogen,dt.FigRSTARLTComp.R33.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)



# FLU-COVID - Long-term competition: COVID 1, FLU 0 + HI - LLI 30-days - R3.3
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP233/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARLTComp30.R33.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARLTComp30.R33.AR$Pathogen<-factor(dt.FigRSTARLTComp30.R33.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp30.R33.AR$IntroductionTime<-factor(dt.FigRSTARLTComp30.R33.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp30.R33.AR<-ggplot(data = dt.FigRSTARLTComp30.R33.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0.65,0.92))
plot(RSTAR.LTComp30.R33.AR)


round(tapply(dt.FigRSTARLTComp30.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp30.R33.AR$Pathogen,dt.FigRSTARLTComp30.R33.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigRSTARLTComp30.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp30.R33.AR$Pathogen,dt.FigRSTARLTComp30.R33.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigRSTARLTComp30.R33.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTComp30.R33.AR$Pathogen,dt.FigRSTARLTComp30.R33.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP233/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOMP33/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R23.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARLTComp30.R33.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARLTComp30.R33.CP$Pathogen<-factor(dt.FigRSTARLTComp30.R33.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTComp30.R33.CP$IntroductionTime<-factor(dt.FigRSTARLTComp30.R33.CP$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTComp30.R33.CP<-ggplot(data = dt.FigRSTARLTComp30.R33.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0.08,0.65))

plot(RSTAR.LTComp30.R33.CP)

tapply(dt.FigRSTARLTComp30.R33.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp30.R33.CP$Pathogen,dt.FigRSTARLTComp30.R33.CP$IntroductionTime),FUN = mean)
tapply(dt.FigRSTARLTComp30.R33.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp30.R33.CP$Pathogen,dt.FigRSTARLTComp30.R33.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigRSTARLTComp30.R33.CP$PeakCases/2500, INDEX = list(dt.FigRSTARLTComp30.R33.CP$Pathogen,dt.FigRSTARLTComp30.R33.CP$IntroductionTime),FUN = quantile, p=0.025 )


RSTAR.COMP.R33.AR <- ggarrange(RSTAR.LTComp.R33.AR + rremove("xlab") + rremove("ylab"), RSTAR.LTComp30.R33.AR + rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                           labels = c("Long-term competition - baseline","Long-term competition - high competion","Long-term competition - baseline","Long-term competition - high competion"),
                           ncol = 2, nrow = 1,
                           common.legend = TRUE, legend = "none",
                           align = "hv", 
                           font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
plot(RSTAR.COMP.R33.AR)
RSTAR.COMP.R33.AR<-annotate_figure(RSTAR.COMP.R33.AR, left = textGrob("Attack Rate", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))


RSTAR.COMP.R33.CP <- ggarrange(RSTAR.LTComp.R33.CP + rremove("xlab") + rremove("ylab"),RSTAR.LTComp30.R33.CP + rremove("xlab") + rremove("ylab"), # remove axis labels from plots
                               labels = c("Long-term competition - baseline","Long-term competition - high competion","Long-term competition - baseline","Long-term competition - high competion"),
                               ncol = 2, nrow = 1,
                               common.legend = TRUE, legend = "none",
                               align = "hv", 
                               font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
plot(RSTAR.COMP.R33.CP)


RSTAR.COMP.R33.CP<-annotate_figure(RSTAR.COMP.R33.CP, left = textGrob("Cases at Peak", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))


ggarrange(RSTAR.COMP.R33.AR,RSTAR.COMP.R33.CP, common.legend = TRUE, ncol = 1, legend = "top")




# FLU-COVID - Long-term co-operation: COVID 1, FLU 0 + HI
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigRSTARLTCoop.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigRSTARLTCoop.AR$Pathogen<-factor(dt.FigRSTARLTCoop.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTCoop.AR$IntroductionTime<-factor(dt.FigRSTARLTCoop.AR$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTCoop.AR<-ggplot(data = dt.FigRSTARLTCoop.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.7))
plot(RSTAR.LTCoop.AR)

tapply(dt.FigRSTARLTCoop.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTCoop.AR$Pathogen,dt.FigRSTARLTCoop.AR$IntroductionTime),FUN = mean)
tapply(dt.FigRSTARLTCoop.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTCoop.AR$Pathogen,dt.FigRSTARLTCoop.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigRSTARLTCoop.AR$FinalSize/2500, INDEX = list(dt.FigRSTARLTCoop.AR$Pathogen,dt.FigRSTARLTCoop.AR$IntroductionTime),FUN = quantile, p=0.025 )


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t235_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t270_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SAMERSTARLCOOP/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R11.3_R21.3_qhqg8.27_t2210_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigRSTARLTCoop.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigRSTARLTCoop.CP$Pathogen<-factor(dt.FigRSTARLTCoop.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigRSTARLTCoop.CP$IntroductionTime<-factor(dt.FigRSTARLTCoop.CP$IntroductionTime, levels = c("0","1","2","3","4"))


RSTAR.LTCoop.CP<-ggplot(data = dt.FigRSTARLTCoop.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.25))
plot(RSTAR.LTCoop.CP)

tapply(dt.FigNHILTCoop.CP$PeakCases/2500, INDEX = list(dt.FigNHILTCoop.CP$Pathogen,dt.FigNHILTCoop.CP$IntroductionTime),FUN = mean)
tapply(dt.FigNHILTCoop.CP$PeakCases/2500, INDEX = list(dt.FigNHILTCoop.CP$Pathogen,dt.FigNHILTCoop.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigNHILTCoop.CP$PeakCases/2500, INDEX = list(dt.FigNHILTCoop.CP$Pathogen,dt.FigNHILTCoop.CP$IntroductionTime),FUN = quantile, p=0.025 )


RSTAR.COOP<- ggarrange(RSTAR.LTCoop.AR + rremove("xlab"), RSTAR.LTCoop.CP + rremove("xlab"), # remove axis labels from plots
                            labels = c("A","B"),
                            ncol = 2, nrow = 1,
                            common.legend = TRUE, legend = "none",
                            align = "hv", 
                            font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))


annotate_figure(RSTAR.COOP, 
                bottom = textGrob("Introduction time influenza", gp = gpar(cex = 1.3)))



#Saved: 950x820
#




# # Different q: ----

# FLU-COVID - Short-term cross-immunity: COVID 1, FLU 0 + HI - q-1
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2Q1/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t248_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t224_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t248_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t2144_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigNHSTCQ1.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHSTCQ1.AR$Pathogen<-factor(dt.FigNHSTCQ1.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHSTCQ1.AR$IntroductionTime<-factor(dt.FigNHSTCQ1.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.STC.Q1.AR<-ggplot(data = dt.FigNHSTCQ1.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.8))+geom_hline(yintercept = 0.1, alpha=0.5,lty=2)
plot(HI.STC.Q1.AR)

round(tapply(dt.FigNHSTCQ1.AR$FinalSize/2500, INDEX = list(dt.FigNHSTCQ1.AR$Pathogen,dt.FigNHSTCQ1.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigNHSTCQ1.AR$FinalSize/2500, INDEX = list(dt.FigNHSTCQ1.AR$Pathogen,dt.FigNHSTCQ1.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigNHSTCQ1.AR$FinalSize/2500, INDEX = list(dt.FigNHSTCQ1.AR$Pathogen,dt.FigNHSTCQ1.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)


temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2Q1/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t248_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t224_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t248_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ1/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t2144_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigNHISTCQ1.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHISTCQ1.CP$Pathogen<-factor(dt.FigNHISTCQ1.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHISTCQ1.CP$IntroductionTime<-factor(dt.FigNHISTCQ1.CP$IntroductionTime, levels = c("0","1","2","3","4"))

HI.STC.Q1.CP<-ggplot(data = dt.FigNHISTCQ1.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.35))
plot(HI.STC.Q1.CP)

# FLU-COVID - Short-term cross-immunity: COVID 1, FLU 0 + HI - q-0.0001
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2Q2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t238_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t219_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t238_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t2114_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))


dt.FigNHSTCQ2.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHSTCQ2.AR$Pathogen<-factor(dt.FigNHSTCQ2.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHSTCQ2.AR$IntroductionTime<-factor(dt.FigNHSTCQ2.AR$IntroductionTime, levels = c("0","1","2","3","4"))


HI.STC.Q2.AR<-ggplot(data = dt.FigNHSTCQ2.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.8))+geom_hline(yintercept = 0.1, alpha=0.5,lty=2)
plot(HI.STC.Q2.AR)

round(tapply(dt.FigNHSTCQ2.AR$FinalSize/2500, INDEX = list(dt.FigNHSTCQ2.AR$Pathogen,dt.FigNHSTCQ2.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigNHSTCQ2.AR$FinalSize/2500, INDEX = list(dt.FigNHSTCQ2.AR$Pathogen,dt.FigNHSTCQ2.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigNHSTCQ2.AR$FinalSize/2500, INDEX = list(dt.FigNHSTCQ2.AR$Pathogen,dt.FigNHSTCQ2.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)



temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2Q2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t238_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t219_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t238_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCIQ2/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t2114_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
}


dt.FigNHISTCQ2.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHISTCQ2.CP$Pathogen<-factor(dt.FigNHISTCQ2.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHISTCQ2.CP$IntroductionTime<-factor(dt.FigNHISTCQ2.CP$IntroductionTime, levels = c("0","1","2","3","4"))

HI.STC.Q2.CP<-ggplot(data = dt.FigNHISTCQ2.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.35))
plot(HI.STC.Q2.CP)

temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.sce<-NULL



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.sce<-c(temp.sce, rep("1",2*length(epi.outbreak)))


dt.FigNHSTC.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigNHSTC.AR$Pathogen<-factor(dt.FigNHSTC.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigNHSTC.AR$IntroductionTime<-factor(dt.FigNHSTC.AR$IntroductionTime, levels = c("0","1","2","3","4"))


round(tapply(dt.FigNHSTC.AR$FinalSize/2500, INDEX = list(dt.FigNHSTC.AR$Pathogen,dt.FigNHSTC.AR$IntroductionTime),FUN = mean),2)
round(tapply(dt.FigNHSTC.AR$FinalSize/2500, INDEX = list(dt.FigNHSTC.AR$Pathogen,dt.FigNHSTC.AR$IntroductionTime),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigNHSTC.AR$FinalSize/2500, INDEX = list(dt.FigNHSTC.AR$Pathogen,dt.FigNHSTC.AR$IntroductionTime),FUN = quantile, p=0.975 ),2)



HI.STC.AR<-ggplot(data = dt.FigNHSTC.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.8))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)
plot(HI.STC.AR)






temp.path.pk<-NULL
temp.t2.pk<-NULL
temp.pk<-NULL
temp.sce.pk<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path2","Path1")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("1",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("2",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("3",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/HOMEISOSCI/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t2150_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("4",2))
  temp.sce.pk<-c(temp.sce.pk, rep("1",2))
}


dt.FigNHISTC.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigNHISTC.CP$Pathogen<-factor(dt.FigNHISTC.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigNHISTC.CP$IntroductionTime<-factor(dt.FigNHISTC.CP$IntroductionTime, levels = c("0","1","2","3","4"))

HI.STC.CP<-ggplot(data = dt.FigNHISTC.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,0.35))
plot(HI.STC.CP)

library(latex2exp)

HOMEISO.STC.Q.AR <- ggarrange(HI.STC.AR + rremove("xlab")+ rremove("ylab"), HI.STC.Q1.AR + rremove("xlab")+ rremove("ylab"), HI.STC.Q2.AR + rremove("ylab")+ rremove("xlab"), # remove axis labels from plots
                           labels = c("Baseline","q_h=q_g","q_h=1e-04"),
                           ncol = 3, nrow = 1,
                           common.legend = TRUE, legend = "top",
#                          align = "hv", 
                            label.x = c(0.1,0.05,0.1),
                          font.label = list(size = 11, color = "black", face = "bold", family = NULL))+ theme(plot.margin = unit(c(-0.2,1,1,1), "lines"))


HOMEISO.STC.Q.AR<-annotate_figure(HOMEISO.STC.Q.AR, left = textGrob("Attack Rate", rot = 90, vjust = 1.5, gp = gpar(cex = 1.3)),
                                  bottom = textGrob("Introduction time influenza", vjust = -0.1, gp = gpar(cex = 1.3)))

plot(HOMEISO.STC.Q.AR)

HOMEISO.STC.Q.CP <- ggarrange(HI.STC.CP + rremove("xlab") + rremove("ylab"), HI.STC.Q1.CP + rremove("xlab") + rremove("ylab"), HI.STC.Q2.CP + rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                              labels = c("Baseline","q_h=q_g","q_h=1e-04"),
                              ncol = 3, nrow = 1,
                           common.legend = TRUE, legend = "none",
                           align = "hv", 
                           label.x = c(0.1,0.05,0.1),
                           font.label = list(size = 11, color = "black", face = "bold", family = NULL, position = "top"))+ theme(plot.margin = unit(c(1.1,1.1,1.1,1.1), "lines"))

HOMEISO.STC.Q.CP<-annotate_figure(HOMEISO.STC.Q.CP, left = textGrob("Cases at Peak", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Introduction time influenza", vjust=-0.25,gp = gpar(cex = 1.3)))

plot(HOMEISO.STC.Q.CP)
HOMEISO.Q<-ggarrange(HOMEISO.STC.Q.AR,HOMEISO.STC.Q.CP, common.legend = TRUE, nrow = 2)
plot(HOMEISO.Q)


# FLU-COVID - GLCTRED - Independent q1
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.cd<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2Q1/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t233_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2Q1/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t233_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2Q1/2/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t233_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t216.5_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t216.5_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t216.5_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t233_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t233_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t233_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t299_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t299_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ1/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t299_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



dt.FigGCRINDQ1.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Reduction"=temp.cd)
dt.FigGCRINDQ1.AR$Pathogen<-factor(dt.FigGCRINDQ1.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigGCRINDQ1.AR$IntroductionTime<-factor(dt.FigGCRINDQ1.AR$IntroductionTime, levels = c("0","1","2","3","4"))
dt.FigGCRINDQ1.AR$Reduction<-factor(dt.FigGCRINDQ1.AR$Reduction, levels = c("2","0","1"))


GCR.IND.Q1.AR<-ggplot(data = dt.FigGCRINDQ1.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,1))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30%","1"="65%")))+geom_hline(yintercept = 0.1, alpha=0.5,lty=2)
plot(GCR.IND.Q1.AR)

round(tapply(dt.FigGCRINDQ1.AR$FinalSize/2500, INDEX = list(dt.FigGCRINDQ1.AR$Pathogen,dt.FigGCRINDQ1.AR$IntroductionTime, dt.FigGCRINDQ1.AR$Reduction),FUN = mean),2)
round(tapply(dt.FigGCRINDQ1.AR$FinalSize/2500, INDEX = list(dt.FigGCRINDQ1.AR$Pathogen,dt.FigGCRINDQ1.AR$IntroductionTime, dt.FigGCRINDQ1.AR$Reduction),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigGCRINDQ1.AR$FinalSize/2500, INDEX = list(dt.FigGCRINDQ1.AR$Pathogen,dt.FigGCRINDQ1.AR$IntroductionTime, dt.FigGCRINDQ1.AR$Reduction),FUN = quantile, p=0.975 ),2)


# FLU-COVID - GLCTRED - Independent q0.0001
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.cd<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2Q2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t227_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2Q2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t227_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDIND2Q2/2/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t227_sigma12_1_sigma21_1_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t213.5_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t213.5_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t213.5_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t227_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t227_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t227_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t281_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t281_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDINDQ2/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t281_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



dt.FigGCRINDQ2.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Reduction"=temp.cd)
dt.FigGCRINDQ2.AR$Pathogen<-factor(dt.FigGCRINDQ2.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigGCRINDQ2.AR$IntroductionTime<-factor(dt.FigGCRINDQ2.AR$IntroductionTime, levels = c("0","1","2","3","4"))
dt.FigGCRINDQ2.AR$Reduction<-factor(dt.FigGCRINDQ2.AR$Reduction, levels = c("2","0","1"))


GCR.IND.Q2.AR<-ggplot(data = dt.FigGCRINDQ2.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,1))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30%","1"="65%")))+geom_hline(yintercept = 0.1, alpha=0.5,lty=2)
plot(GCR.IND.Q2.AR)

round(tapply(dt.FigGCRINDQ2.AR$FinalSize/2500, INDEX = list(dt.FigGCRINDQ2.AR$Pathogen,dt.FigGCRINDQ2.AR$IntroductionTime, dt.FigGCRINDQ2.AR$Reduction),FUN = mean),2)
round(tapply(dt.FigGCRINDQ2.AR$FinalSize/2500, INDEX = list(dt.FigGCRINDQ2.AR$Pathogen,dt.FigGCRINDQ2.AR$IntroductionTime, dt.FigGCRINDQ2.AR$Reduction),FUN = quantile, p=0.025 ),2)
round(tapply(dt.FigGCRINDQ2.AR$FinalSize/2500, INDEX = list(dt.FigGCRINDQ2.AR$Pathogen,dt.FigGCRINDQ2.AR$IntroductionTime, dt.FigGCRINDQ2.AR$Reduction),FUN = quantile, p=0.975 ),2)



GCR.IND.AR<-ggplot(data = dt.FigGCRIND.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,1))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30% reduction","1"="65% reduction")))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)

ggarrange(GCR.IND.Q2.AR,GCR.IND.Q1.AR,GCR.IND.AR, ncol = 3, labels = c("A","B","C"), common.legend = TRUE, legend = "top", label.x = 0.05)


# FLU-COVID - Cross-immunity during infection  - q1
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.cd<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2Q1/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t233_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2Q1/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t233_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2Q1/2/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1_t233_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t216.5_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t216.5_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t216.5_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t233_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t233_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t233_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t299_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t299_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ1/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t299_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



dt.FigGCRSCIQ1.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Reduction"=temp.cd)
dt.FigGCRSCIQ1.AR$Pathogen<-factor(dt.FigGCRSCIQ1.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigGCRSCIQ1.AR$IntroductionTime<-factor(dt.FigGCRSCIQ1.AR$IntroductionTime, levels = c("0","1","2","3","4"))
dt.FigGCRSCIQ1.AR$Reduction<-factor(dt.FigGCRSCIQ1.AR$Reduction, levels = c("2","0","1"))


GCR.SCI.Q1.AR<-ggplot(data = dt.FigGCRSCIQ1.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,1))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30%","1"="65%")))+geom_hline(yintercept = 0.1, alpha=0.5,lty=2)
plot(GCR.SCI.Q1.AR)

# FLU-COVID - Cross-immunity during infection  - q 1e-04
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL
temp.cd<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2Q2/0/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t227_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2Q2/1/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t227_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCI2Q2/2/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg1e-04_t227_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path2","Path1"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("1",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t213.5_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t213.5_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t213.5_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("2",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t227_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t227_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t227_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("3",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t281_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("0",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t281_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("1",2*length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/GLCTREDSCIQ2/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1e-04_t281_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("4",2*length(epi.outbreak)))
temp.cd<-c(temp.cd, rep("2",2*length(epi.outbreak)))



dt.FigGCRSCIQ2.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Reduction"=temp.cd)
dt.FigGCRSCIQ2.AR$Pathogen<-factor(dt.FigGCRSCIQ2.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigGCRSCIQ2.AR$IntroductionTime<-factor(dt.FigGCRSCIQ2.AR$IntroductionTime, levels = c("0","1","2","3","4"))
dt.FigGCRSCIQ2.AR$Reduction<-factor(dt.FigGCRSCIQ2.AR$Reduction, levels = c("2","0","1"))


GCR.SCI.Q2.AR<-ggplot(data = dt.FigGCRSCIQ2.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,1))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30%","1"="65%")))+geom_hline(yintercept = 0.1, alpha=0.5,lty=2)
plot(GCR.SCI.Q2.AR)

GCR.SCI.AR<-ggplot(data = dt.FigGCRSCI.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_boxplot()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Precendent","Concurrent","Ascending","Peak","Subsequent"))+labs(fill="Disease")+ylim(c(0,1))+facet_grid(~Reduction,  scales = "fixed", labeller = labeller(Reduction=c("2"="No reduction" ,"0"="30% reduction","1"="65% reduction")))+geom_hline(yintercept = 0.1,lty=2,alpha=0.5)

ggarrange(GCR.SCI.Q2.AR, GCR.SCI.Q1.AR, GCR.SCI.AR, ncol = 3, labels = c("A","B","C"), common.legend = TRUE, label.x = 0.01)

#



######################################################
# Appendix
#####################################################
# HOUSEHOLD SIZE DENSITY ----
freq<-c(rep(1,163186),rep(2,151288),rep(3,72304),rep(4,60656),rep(5,22343),rep(6,6889),rep(7,2097))
hhsz.dens<-data.frame("Freq"=freq)

ggplot(hhsz.dens, aes(hhsz.dens$Freq))+geom_histogram(aes(y=..density..),binwidth = 1, color="#e9ecef", fill="steelblue")+
  theme_minimal()+theme(
    legend.key.size = unit(0.8, "cm"),
    legend.key.width = unit(0.8,"cm"),
    legend.title = element_text(size=17, face = "bold"),
    legend.text = element_text(size=16),
    axis.text.x = element_text(size=16),
    axis.title = element_text(size=20),
    axis.text.y = element_text(size=16),
    text = element_text(size=20, face = "bold"),
  ) + xlab("Household Size")+ylab("Density")+scale_x_continuous(breaks = 1:7)




# Fig X Independent pathogen and no compliance ----
temp.path<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t254_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigXa.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigXa.AR$Pathogen<-factor(dt.FigXa.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigXa.AR$IntroductionTime<-factor(dt.FigXa.AR$IntroductionTime, levels = c("0","13","26","39"))


aa101<-ggplot(data = dt.FigXa.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,1))

temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(123)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t218_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t236_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPIND/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t254_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc20.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.FigXa.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigXa.CP$Pathogen<-factor(dt.FigXa.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigXa.CP$IntroductionTime<-factor(dt.FigXa.CP$IntroductionTime, levels = c("0","13","26","39"))


aa31<-ggplot(data = dt.FigXa.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","FLU"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ylim(c(0,0.5))

ggarrange(aa101,aa31, nrow = 1, ncol = 2, labels = c("A","B"), common.legend = TRUE)

print("A")
tapply(dt.FigXa.AR$FinalSize /2500, INDEX = list(dt.FigXa.AR$Pathogen,dt.FigXa.AR$IntroductionTime),FUN = mean)
tapply(dt.FigXa.AR$FinalSize /2500, INDEX = list(dt.FigXa.AR$Pathogen,dt.FigXa.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigXa.AR$FinalSize/2500, INDEX = list(dt.FigXa.AR$Pathogen,dt.FigXa.AR$IntroductionTime),FUN = quantile, p=0.025 )

print("B")
tapply(dt.FigXa.CP$PeakCases /2500, INDEX = list(dt.FigXa.CP$Pathogen,dt.FigXa.CP$IntroductionTime),FUN = mean)
tapply(dt.FigXa.CP$PeakCases/2500, INDEX = list(dt.FigXa.CP$Pathogen,dt.FigXa.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigXa.CP$PeakCases/2500, INDEX = list(dt.FigXa.CP$Pathogen,dt.FigXa.CP$IntroductionTime),FUN = quantile, p=0.025 )


# FIG S1 Pathogens and their viral characteristics ----
# Pathogens and their viral characteristics
time.points1<-seq(0,12,0.1)
time.points2<-seq(0,8,0.1)
inf.covid19<-dgamma(time.points1,shape = 12, rate = 2.08)/ (pgamma(15,shape = 2,rate = 2.08))
inf.fluA<-1.001592 *dgamma(time.points2,shape = 4.604016,scale = 0.5922317)

Inf.Measure<-data.frame(x=c(time.points1,time.points2),y=c(inf.covid19,inf.fluA) , Disease=c(rep("COVID-19",length(time.points1)),rep("FLU-A",length(time.points2))), Quantity=rep("InfectiousnessMeasure",length(time.points1)+length(time.points2)))

dt<-data.frame(x=Inf.Measure$x, y=Inf.Measure$y,Disease=c(as.character(Inf.Measure$Disease)), Quantity=c(as.character(Inf.Measure$Quantity)))
pt<-data.frame(x=c(0,0), y=c(0,0), Quantity=c("InfectiousnessMeasure"))

ggplot(dt, aes(x=x, y=y))+geom_line(data = Inf.Measure, aes(x=x, y=y, col=Disease),lwd=1.3,lty=1)+
  theme_classic()+theme(
    legend.key.size = unit(0.8, "cm"),
    legend.key.width = unit(0.8,"cm"),
    legend.title = element_text(size=14, face = "bold"),
    legend.text = element_text(size=13),
    #    legend.position = "top",
    axis.text.x = element_text(size=13),
    axis.title = element_text(size=14),
    axis.text.y = element_text(size=13),
    text = element_text(size=17, face = "bold"),
  ) + xlab("Days since infection")+ylab("Infectiousness Measure")+scale_color_manual(values=c("#56B4E9","#E69F00"),labels=c("COVID-19","Influenza"))+geom_vline(xintercept = 11/2.08,colour="#56B4E9",lty=2,lwd=0.8)+geom_vline(xintercept = 2,colour="#E69F00",lty=2,lwd=0.8)+ylim(c(0,0.4))+xlim(c(0,12.5))+labs(color="Disease:")

#

#saved 650*350



# Fig S2 Decrease in local contact rate ----
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISo/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS2<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS2$Pathogen<-factor(dt.FigS2$Pathogen, levels = c("Path1","Path2"))
dt.FigS2$IntroductionTime<-factor(dt.FigS2$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS2$Compliance<-factor(dt.FigS2$Compliance, levels = c("1","0.5","0"))


aa101<-ggplot(data = dt.FigS2, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.7))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("0"="Complete Isolation","0.5"="50% decrease","1"="No decrease")))

plot(aa101)

#Figure 1
print("A")
tapply(dt.FigS2$FinalSize/2500, INDEX = list(dt.FigS2$Pathogen,dt.FigS2$IntroductionTime, dt.FigS2$Compliance ),FUN = mean)
tapply(dt.FigS2$FinalSize/2500, INDEX = list(dt.FigS2$Pathogen,dt.FigS2$IntroductionTime, dt.FigS2$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS2$FinalSize/2500, INDEX = list(dt.FigS2$Pathogen,dt.FigS2$IntroductionTime, dt.FigS2$Compliance),FUN = quantile, p=0.025 )




# Fig S3 Competition profiles and length ----
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/13/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/16/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/19/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/22/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/25/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/28/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/31/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/34/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS3a<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS3a$Pathogen<-factor(dt.FigS3a$Pathogen, levels = c("Path1","Path2"))
dt.FigS3a$IntroductionTime<-factor(dt.FigS3a$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS3a$Compliance<-factor(dt.FigS3a$Compliance, levels = c("1","3","2"))


aa101<-ggplot(data = dt.FigS3a, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.69))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="Constant","2"="Decreasing2","3"="Decreasing1")))
plot(aa101)

temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/18/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/19/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/20/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/27/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/28/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/29/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS3b<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS3b$Pathogen<-factor(dt.FigS3b$Pathogen, levels = c("Path1","Path2"))
dt.FigS3b$IntroductionTime<-factor(dt.FigS3b$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS3b$Compliance<-factor(dt.FigS3b$Compliance, levels = c("1","2","3"))


aa102<-ggplot(data = dt.FigS3b, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.69))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="No reduced susceptibility","2"="10 days reduced susceptibility","3"="30 days reduced susceptibility")))
plot(aa102)

ggarrange(aa101,aa102, nrow = 2, ncol = 1, labels = c("A","B"), common.legend = TRUE)

#Saved: 950 x 820
#Table
#Figure 1
print("A")
tapply(dt.FigS3a$FinalSize/2500, INDEX = list(dt.FigS3a$Pathogen,dt.FigS3a$IntroductionTime, dt.FigS3a$Compliance ),FUN = mean)
tapply(dt.FigS3a$FinalSize/2500, INDEX = list(dt.FigS3a$Pathogen,dt.FigS3a$IntroductionTime, dt.FigS3a$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS3a$FinalSize/2500, INDEX = list(dt.FigS3a$Pathogen,dt.FigS3a$IntroductionTime, dt.FigS3a$Compliance),FUN = quantile, p=0.025 )
#Figure 2
print("B")
tapply(dt.FigS3b$FinalSize/2500, INDEX = list(dt.FigS3b$Pathogen,dt.FigS3b$IntroductionTime, dt.FigS3b$Compliance),FUN = mean)
tapply(dt.FigS3b$FinalSize/2500, INDEX = list(dt.FigS3b$Pathogen,dt.FigS3b$IntroductionTime, dt.FigS3b$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS3b$FinalSize/2500, INDEX = list(dt.FigS3b$Pathogen,dt.FigS3b$IntroductionTime, dt.FigS3b$Compliance),FUN = quantile, p=0.025 )

# Fig S4 Cooperation profiles and length ----
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/13/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/16/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/22/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/19/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/25/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/31/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/28/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/34/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS4a<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS4a$Pathogen<-factor(dt.FigS4a$Pathogen, levels = c("Path1","Path2"))
dt.FigS4a$IntroductionTime<-factor(dt.FigS4a$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS4a$Compliance<-factor(dt.FigS4a$Compliance, levels = c("1","3","2"))


aa101<-ggplot(data = dt.FigS4a, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.8))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="Constant","2"="Decreasing2","3"="Decreasing1")))
plot(aa101)

temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/12/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/13/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/14/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/21/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/22/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/23/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/30/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/31/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/32/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("3",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS4b<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS4b$Pathogen<-factor(dt.FigS4b$Pathogen, levels = c("Path1","Path2"))
dt.FigS4b$IntroductionTime<-factor(dt.FigS4b$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS4b$Compliance<-factor(dt.FigS4b$Compliance, levels = c("1","2","3"))


aa102<-ggplot(data = dt.FigS4b, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c( "#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.8))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="No increased susceptibility ","2"="10 days increased susceptibility","3"="30 days increased susceptibility")))
plot(aa102)

ggarrange(aa101,aa102, nrow = 2, ncol = 1, labels = c("A","B"), common.legend = TRUE)

#Saved: 950 x 820
#Table
#Figure 1
print("A")
tapply(dt.FigS4a$FinalSize/2500, INDEX = list(dt.FigS4a$Pathogen,dt.FigS4a$IntroductionTime, dt.FigS4a$Compliance ),FUN = mean)
tapply(dt.FigS4a$FinalSize/2500, INDEX = list(dt.FigS4a$Pathogen,dt.FigS4a$IntroductionTime, dt.FigS4a$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS4a$FinalSize/2500, INDEX = list(dt.FigS4a$Pathogen,dt.FigS4a$IntroductionTime, dt.FigS4a$Compliance),FUN = quantile, p=0.025 )
#Figure 2
print("B")
tapply(dt.FigS4b$FinalSize/2500, INDEX = list(dt.FigS4b$Pathogen,dt.FigS4b$IntroductionTime, dt.FigS4b$Compliance),FUN = mean)
tapply(dt.FigS4b$FinalSize/2500, INDEX = list(dt.FigS4b$Pathogen,dt.FigS4b$IntroductionTime, dt.FigS4b$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS4b$FinalSize/2500, INDEX = list(dt.FigS4b$Pathogen,dt.FigS4b$IntroductionTime, dt.FigS4b$Compliance),FUN = quantile, p=0.025 )





# Fig S5 Value of Q ----
#Cross-protection co-infection
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/18/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/27/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS5a<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS5a$Pathogen<-factor(dt.FigS5a$Pathogen, levels = c("Path1","Path2"))
dt.FigS5a$IntroductionTime<-factor(dt.FigS5a$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS5a$Compliance<-factor(dt.FigS5a$Compliance, levels = c("1","2"))


aa102<-ggplot(data = dt.FigS5a, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.69))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="qh=8.27qg","2"="qh=qg")))
plot(aa102)

print("A")
tapply(dt.FigS5a$FinalSize/2500, INDEX = list(dt.FigS5a$Pathogen,dt.FigS5a$IntroductionTime, dt.FigS5a$Compliance ),FUN = mean)
tapply(dt.FigS5a$FinalSize/2500, INDEX = list(dt.FigS5a$Pathogen,dt.FigS5a$IntroductionTime, dt.FigS5a$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS5a$FinalSize/2500, INDEX = list(dt.FigS5a$Pathogen,dt.FigS5a$IntroductionTime, dt.FigS5a$Compliance),FUN = quantile, p=0.025 )


#Competition - 10 days
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/19/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHIComp/28/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICompQ/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS5b<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS5b$Pathogen<-factor(dt.FigS5b$Pathogen, levels = c("Path1","Path2"))
dt.FigS5b$IntroductionTime<-factor(dt.FigS5b$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS5b$Compliance<-factor(dt.FigS5b$Compliance, levels = c("1","2"))


ab102<-ggplot(data = dt.FigS5b, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.69))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="qh=8.27qg","2"="qh=qg")))
plot(ab102)

print("B")
tapply(dt.FigS5b$FinalSize/2500, INDEX = list(dt.FigS5b$Pathogen,dt.FigS5b$IntroductionTime, dt.FigS5b$Compliance ),FUN = mean)
tapply(dt.FigS5b$FinalSize/2500, INDEX = list(dt.FigS5b$Pathogen,dt.FigS5b$IntroductionTime, dt.FigS5b$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS5b$FinalSize/2500, INDEX = list(dt.FigS5b$Pathogen,dt.FigS5b$IntroductionTime, dt.FigS5b$Compliance),FUN = quantile, p=0.025 )



#Cooperation - 10 days - 
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))



load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOGQ/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/13/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOGQ/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t225_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/22/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOGQ/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t250_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOG/31/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg8.27_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopOGQ/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t275_sigma12_0_sigma21_0_alpha10.5_alpha20.33_rho10.69_rho20.67_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("2",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS5c<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS5c$Pathogen<-factor(dt.FigS5c$Pathogen, levels = c("Path1","Path2"))
dt.FigS5c$IntroductionTime<-factor(dt.FigS5c$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS5c$Compliance<-factor(dt.FigS5c$Compliance, levels = c("1","2"))


ac102<-ggplot(data = dt.FigS5c, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.7))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("1"="qh=8.27qg","2"="qh=qg")))
plot(ac102)

print("C")
tapply(dt.FigS5c$FinalSize/2500, INDEX = list(dt.FigS5c$Pathogen,dt.FigS5c$IntroductionTime, dt.FigS5c$Compliance ),FUN = mean)
tapply(dt.FigS5c$FinalSize/2500, INDEX = list(dt.FigS5c$Pathogen,dt.FigS5c$IntroductionTime, dt.FigS5c$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS5c$FinalSize/2500, INDEX = list(dt.FigS5c$Pathogen,dt.FigS5c$IntroductionTime, dt.FigS5c$Compliance),FUN = quantile, p=0.025 )


# Home iso - ctc rate
temp.path<-NULL
temp.hcrc<-NULL
temp.t2<-NULL
temp.fs<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/0/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/1/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/2/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t20_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/3/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/4/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/5/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t225_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("25",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/6/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/7/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/8/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t250_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("50",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/9/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/10/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed0.5_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("0.5",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHomISoQ/11/MP_Synth_nVertex2500_nNetw100COVID-19_&_FLU-A_R13.3_R21.3_qhqg1_t275_sigma12_1_sigma21_1_alpha10.5_alpha20.33_rho10.69_rho20.67_lli11_lli21_NetSynth_CtcRed1_PImm0_tSeed1000_bc11_bc20.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.hcrc<-c(temp.hcrc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("75",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS5d<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.hcrc)
dt.FigS5d$Pathogen<-factor(dt.FigS5d$Pathogen, levels = c("Path1","Path2"))
dt.FigS5d$IntroductionTime<-factor(dt.FigS5d$IntroductionTime, levels = c("0","25","50","75"))
dt.FigS5d$Compliance<-factor(dt.FigS5d$Compliance, levels = c("1","0.5","0"))


aa101<-ggplot(data = dt.FigS5d, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#56B4E9","#E69F00"), labels=c("COVID-19","Influenza"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time influenza")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.7))+facet_grid(~Compliance,  scales = "fixed", labeller = labeller(Compliance=c("0"="Complete Isolation","0.5"="50% decrease","1"="No decrease")))

plot(aa101)

print("A")
tapply(dt.FigS5d$FinalSize/2500, INDEX = list(dt.FigS5d$Pathogen,dt.FigS5d$IntroductionTime, dt.FigS5d$Compliance ),FUN = mean)
tapply(dt.FigS5d$FinalSize/2500, INDEX = list(dt.FigS5d$Pathogen,dt.FigS5d$IntroductionTime, dt.FigS5d$Compliance),FUN = quantile, p=0.975 )
tapply(dt.FigS5d$FinalSize/2500, INDEX = list(dt.FigS5d$Pathogen,dt.FigS5d$IntroductionTime, dt.FigS5d$Compliance),FUN = quantile, p=0.025 )




# FIG S6 Long-term Competition + Cooperation (Influenza seeded during COVID-19 epidemic) ----
temp.path<-NULL
temp.bc<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/8/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/17/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t212.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/26/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/35/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t237.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS6a.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2)
dt.FigS6a.AR$Pathogen<-factor(dt.FigS6a.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigS6a.AR$IntroductionTime<-factor(dt.FigS6a.AR$IntroductionTime, levels = c("0","13","26","39"))


aa101<-ggplot(data = dt.FigS6a.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#E69F00","#56B4E9"), labels=c("Influenza","COVID-19"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time COVID-19")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.9))


# FLU-COVID - Prevalence curves over time
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(12345)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/8/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/17/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t212.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/26/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoopSM/35/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t237.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli10.5_lli20.5_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path1","Path2"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","13","26","39"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","13","26","39"),z=c(0,25,50,75))

aa21<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","FLU"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "grey",
                                  size = 0.2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=10),
  legend.position = "top",
  legend.title = element_text(size=11, face = "bold"),
  axis.text.x = element_text(size=10, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=11),
  axis.text.y = element_text(size=10),
  strip.text = element_text(size = 11)
)+ylab(("Prevalence"))+xlab("Introduction time COVID-19")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Beginning","13"="Ascending","26"="Peak","39"="Descending")))+ylim(c(0,0.3))

dt.FigS6a.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigS6a.CP$Pathogen<-factor(dt.FigS6a.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigS6a.CP$IntroductionTime<-factor(dt.FigS6a.CP$IntroductionTime, levels = c("0","13","26","39"))


aa31<-ggplot(data = dt.FigS6a.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#E69F00","#56B4E9"),labels=c("Influenza","COVID-19"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time COVID-19")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ ylim(c(0,0.25))

# FLU-COVID -  2 co-op 

temp.path<-NULL
temp.bc<-NULL
temp.t2<-NULL
temp.fs<-NULL


load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/2/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("0",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("0",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/8/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t212.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("13",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/14/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("26",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/20/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t237.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:length(epi.outbreak)){
  temp.fs<-c(temp.fs,epi.outbreak[[i]]$FinalSize$FinalSize1,epi.outbreak[[i]]$FinalSize$FinalSize2)
}
temp.bc<-c(temp.bc,rep("1",2*length(epi.outbreak)))
temp.t2<-c(temp.t2, rep("39",2*length(epi.outbreak)))
temp.path<-c(temp.path,rep(c("Path1","Path2"),length(epi.outbreak)))

dt.FigS6b.AR<-data.frame("FinalSize"=temp.fs,"Pathogen"=temp.path,"IntroductionTime"=temp.t2, "Compliance"=temp.bc)
dt.FigS6b.AR$Pathogen<-factor(dt.FigS6b.AR$Pathogen, levels = c("Path1","Path2"))
dt.FigS6b.AR$IntroductionTime<-factor(dt.FigS6b.AR$IntroductionTime, levels = c("0","13","26","39"))


aa102<-ggplot(data = dt.FigS6b.AR, aes(x=IntroductionTime, fill=Pathogen,y=FinalSize/2500))+geom_violin()+scale_fill_manual(values = c("#E69F00","#56B4E9"),labels=c("Influenza","COVID-19"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Introduction time COVID-19")+scale_x_discrete(labels=c("Beginning","Ascending ","Peak","Descending"))+labs(fill="Disease")+ylim(c(0,0.9))


# FLU-COVID - Prevalence curves over time
temp.path<-NULL
temp.path.pk<-NULL
temp.t2<-NULL
temp.t2.pk<-NULL
temp.prev<-NULL
temp.time<-NULL
temp.nc<-NULL
ncurves<-10
a<-1:(8*ncurves)
a<-as.character(a)
set.seed(123)
temp.pk<-NULL

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/2/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t20_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("0",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("0",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/8/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t212.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("13",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("13",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/14/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t225_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("26",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("26",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

load("~/Library/CloudStorage/GoogleDrive-andrea.torneri@uhasselt.be/My Drive/Work/PhD/Co-infection/Sims/SIMPAPHICoop/20/MP_Synth_nVertex2500_nNetw100FLU-A_&_COVID-19_R11.3_R23.3_qhqg8.27_t237.5_sigma12_0_sigma21_0_alpha10.33_alpha20.5_rho10.67_rho20.69_lli12_lli22_NetSynth_CtcRed1_PImm0_tSeed1000_bc10_bc21.RData")
for (i in 1:ncurves){
  epi.det<-epi.outbreak[[sample(1:length(epi.outbreak),1)]]$epi.details
  temp.time<-c(temp.time,epi.det$Days,epi.det$Days)
  temp.prev<-c(temp.prev,epi.det$Prevalence1,epi.det$Prevalence2)
  temp.t2<-c(temp.t2, rep("39",2*length(epi.det$Days)))
  temp.path<-c(temp.path,c(rep("Path1",length(epi.det$Days)),rep("Path2",length(epi.det$Days))))
  temp.nc<-c(temp.nc, rep(a[i],2*length(epi.det$Days)))
}
#temp.pk<-NULL
for (i in 1:length(epi.outbreak)){
  temp.pk<-c(temp.pk, epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence1,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence2)
  temp.path.pk<-c(temp.path.pk,rep(c("Path1","Path2")))
  temp.t2.pk<-c(temp.t2.pk, rep("39",2))
}
quantile(temp.pk,0.025)
quantile(temp.pk,0.975)

dt.Prev<-data.frame("Prevalence"=temp.prev,"Days"=temp.time,"IntroductionTime"=temp.t2, "Pathogen"=temp.path, "TempC"=temp.nc)
dt.Prev$Pathogen<-factor(dt.Prev$Pathogen, levels = c("Path2","Path1"))
dt.Prev$IntroductionTime<-factor(dt.Prev$IntroductionTime, levels = c("0","13","26","39"))

dt.1<-dt.Prev[which(dt.Prev$TempC==1),]
dt.2<-dt.Prev[which(dt.Prev$TempC==2),]
dt.3<-dt.Prev[which(dt.Prev$TempC==3),]
dt.4<-dt.Prev[which(dt.Prev$TempC==4),]
dt.5<-dt.Prev[which(dt.Prev$TempC==5),]
dt.6<-dt.Prev[which(dt.Prev$TempC==6),]
dt.7<-dt.Prev[which(dt.Prev$TempC==7),]
dt.8<-dt.Prev[which(dt.Prev$TempC==8),]
dt.9<-dt.Prev[which(dt.Prev$TempC==9),]
dt.10<-dt.Prev[which(dt.Prev$TempC==10),]

library(dplyr)
vl.dt<-data.frame(IntroductionTime=c("0","13","26","39"),z=c(0,12.5,25,37.5))

aa22<-ggplot(data = dt.Prev, aes(x=Days,y=Prevalence/2500, col=Pathogen))+geom_line(data=dt.1)+geom_line(data = dt.2)+geom_line(data = dt.3)+geom_line(data = dt.4)+geom_line(data = dt.5)+geom_line(data = dt.6)+geom_line(data = dt.7)+geom_line(data = dt.8)+geom_line(data = dt.9)+geom_line(data = dt.10)+scale_color_manual(values = c("#56B4E9","#E69F00"),labels=c("COVID-19","FLU"))+geom_vline(aes(xintercept=z),vl.dt, lty=2)+ theme(
  panel.background = element_rect(fill = "white",
                                  colour = "grey",
                                  size = 0.2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Prevalence"))+xlab("Introduction time COVID-19")+labs(colour="Disease")+facet_grid(~IntroductionTime,  scales = "fixed", labeller = labeller(IntroductionTime=c("0"="Beginning","13"="Ascending","26"="Peak","39"="Descending")))+ylim(c(0,0.3))

dt.FigS6b.CP<-data.frame("PeakCases"=temp.pk,"IntroductionTime"=temp.t2.pk, "Pathogen"=temp.path.pk)
dt.FigS6b.CP$Pathogen<-factor(dt.FigS6b.CP$Pathogen, levels = c("Path1","Path2"))
dt.FigS6b.CP$IntroductionTime<-factor(dt.FigS6b.CP$IntroductionTime, levels = c("0","13","26","39"))


aa32<-ggplot(data = dt.FigS6b.CP, aes(x=IntroductionTime,y=PeakCases/2500, fill=Pathogen))+geom_violin()+scale_fill_manual(values = c("#E69F00","#56B4E9"),labels=c("Influenza","COVID-19"))+stat_summary(fun.y=mean, geom="point", shape=15, size=2, color="black", fill="black", aes(group=Pathogen),position=position_dodge(.9))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Cases at Peak"))+xlab("Introduction time COVID-19")+scale_x_discrete(labels=c("Beginning","Ascending " ,"Peak", "Descending"))+labs(fill="Disease")+ylim(c(0,0.25))

ggarrange(aa101,aa102,aa31,aa32, nrow = 2, ncol = 2, labels = c("A","B","C","D"), common.legend = TRUE)

#Saved: 950 x 820
#Table
#Figure 1
print("A")
tapply(dt.FigS6a.AR$FinalSize/2500, INDEX = list(dt.FigS6a.AR$Pathogen,dt.FigS6a.AR$IntroductionTime),FUN = mean)
tapply(dt.FigS6a.AR$FinalSize/2500, INDEX = list(dt.FigS6a.AR$Pathogen,dt.FigS6a.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigS6a.AR$FinalSize/2500, INDEX = list(dt.FigS6a.AR$Pathogen,dt.FigS6a.AR$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 2
print("B")
tapply(dt.FigS6b.AR$FinalSize/2500, INDEX = list(dt.FigS6b.AR$Pathogen,dt.FigS6b.AR$IntroductionTime),FUN = mean)
tapply(dt.FigS6b.AR$FinalSize/2500, INDEX = list(dt.FigS6b.AR$Pathogen,dt.FigS6b.AR$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigS6b.AR$FinalSize/2500, INDEX = list(dt.FigS6b.AR$Pathogen,dt.FigS6b.AR$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 4
print("C")
tapply(dt.FigS6a.CP$PeakCases/2500, INDEX = list(dt.FigS6a.CP$Pathogen,dt.FigS6a.CP$IntroductionTime),FUN = mean)
tapply(dt.FigS6a.CP$PeakCases/2500, INDEX = list(dt.FigS6a.CP$Pathogen,dt.FigS6a.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigS6a.CP$PeakCases/2500, INDEX = list(dt.FigS6a.CP$Pathogen,dt.FigS6a.CP$IntroductionTime),FUN = quantile, p=0.025 )
#Figure 4
print("D")
tapply(dt.FigS6b.CP$PeakCases/2500, INDEX = list(dt.FigS6b.CP$Pathogen,dt.FigS6b.CP$IntroductionTime),FUN = mean)
tapply(dt.FigS6b.CP$PeakCases/2500, INDEX = list(dt.FigS6b.CP$Pathogen,dt.FigS6b.CP$IntroductionTime),FUN = quantile, p=0.975 )
tapply(dt.FigS6b.CP$PeakCases/2500, INDEX = list(dt.FigS6b.CP$Pathogen,dt.FigS6b.CP$IntroductionTime),FUN = quantile, p=0.025 )







#Extinction Fig4
1-(length(which(dt.Fig4a.AR$FinalSize[which(dt.Fig4a.AR$Pathogen=="Path2" & dt.Fig4a.AR$IntroductionTime=="0")]/2500<0.1)))/100
1-(length(which(dt.Fig4a.AR$FinalSize[which(dt.Fig4a.AR$Pathogen=="Path2" & dt.Fig4a.AR$IntroductionTime=="13")]/2500<0.1)))/100
1-(length(which(dt.Fig4a.AR$FinalSize[which(dt.Fig4a.AR$Pathogen=="Path2" & dt.Fig4a.AR$IntroductionTime=="26")]/2500<0.1)))/100
1-(length(which(dt.Fig4a.AR$FinalSize[which(dt.Fig4a.AR$Pathogen=="Path2" & dt.Fig4a.AR$IntroductionTime=="39")]/2500<0.1)))/100

#Extinction S6
1-(length(which(dt.FigS6a.AR$FinalSize[which(dt.FigS6a.AR$Pathogen=="Path1" & dt.FigS6a.AR$IntroductionTime=="0")]/2500<0.1)))/100
1-(length(which(dt.FigS6a.AR$FinalSize[which(dt.FigS6a.AR$Pathogen=="Path1" & dt.FigS6a.AR$IntroductionTime=="13")]/2500<0.1)))/100
1-(length(which(dt.FigS6a.AR$FinalSize[which(dt.FigS6a.AR$Pathogen=="Path1" & dt.FigS6a.AR$IntroductionTime=="26")]/2500<0.1)))/100
1-(length(which(dt.FigS6a.AR$FinalSize[which(dt.FigS6a.AR$Pathogen=="Path1" & dt.FigS6a.AR$IntroductionTime=="39")]/2500<0.1)))/100
