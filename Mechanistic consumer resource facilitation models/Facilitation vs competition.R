### continue on from Facilitation modification simulation script

facil_comp_isolation <- function(t,Compos ,Dose , r,a, pars1,Y) {
  parstt <- pars1
  with( as.list( c(pars1) ), {
    parstt[names(parstt)=="K_N"] <- (a)*pars1[names(parstt)=="K_N"] 
    parstt[names(parstt)=="K_A"] <- pars1[names(parstt)=="K_A"]
    dd_t <- Y[time==t&Compostition==Compos &DoseNum==Dose][rho==r][alpha==a]
    dd_t[State=="N"]$y0
    N <- dd_t[State=="N"]$y0
    A <- dd_t[State=="A"]$y0
    Z_N<- dd_t[State=="Z_N"]$y0
    Z_A<- dd_t[State=="Z_A"]$y0
    S_N<- dd_t[State=="S_N"]$y0
    S_A<- dd_t[State=="S_A"]$y0
    TotA<- dd_t[State=="TotA"]$y0
    TotN<- dd_t[State=="TotN"]$y0
    
    E <- dd_t[State=="E"]$y0
    COMPA <-  (A + Z_A + S_A )/K_A # Resource limitation
    COMPN <-  (N + Z_N + S_N )/K_N 
    
    COMPtot <- 1/(1- (COMPA +COMPN ))
    FACIL <-  1+(E/(1+c_N*E))
    CF <- COMPtot/FACIL

   # COMPtot <- (1- (COMPA +COMPN ))
   # FACIL <-  1+(E/(1+c_N*E))
   # CFN <- r_N*COMPtot*FACIL
   # CFA <- r_A*COMPtot*FACIL
    data.table(time=t,
               Compostition=Compos,
               DoseNum=Dose,rho=r,alpha=a,
               COMPA=COMPA, COMPN=COMPN,COMPtot=COMPtot, FACIL=FACIL,CF=CF,r=r_N,TotN=TotN,TotA=TotA)#CFN=CFN,CFA=CFA)
  } )
}


facil_comp_isolation(t=1,Compos="polyculture" ,Dose=0 , r=1,a=1, pars1=pars1, Y=res)


key1<- expand.grid(time=1:21,DoseNum=res$DoseNum%>%unique)
CF<-rbindlist( lapply(1:nrow(key1) ,function(i){
  tm<-key1[i,]$time
  ds<-key1[i,]$DoseNum
  y<-facil_comp_isolation(t=tm,Compos="polyculture" ,Dose=ds , r=1,a=1, pars1=pars1, Y=res)
  y$whichSpp<- c("both")
  return(y)
}))
CF2<-rbindlist( lapply(1:nrow(key1) ,function(i){
  tm<-key1[i,]$time
  ds<-key1[i,]$DoseNum
  y<-facil_comp_isolation(t=tm,Compos="monoculture" ,Dose=ds , r=1,a=1, pars1=pars1, Y=res)
  y$whichSpp<- c("Nmono","Amono")
  return(y)
}))

CFmono_co<-rbind(CF,CF2)[]

CFlong<-data.table(CFmono_co%>%gather(Process,Value,COMPA:FACIL))
# ggplot(CFlong,aes(y=Value,x=time,col=Process))+geom_line(size=4)
# ggplot(CF,aes(y=COMPtot,x=time))+geom_point()
# ggplot(CF,aes(y=CF,x=time,col=DoseNum))+geom_point()+facet_wrap(~DoseNum)
# 
# ggplot(CF,aes(y=COMPtot,x=FACIL,col=log(1+DoseNum),group=DoseNum))+geom_path()#+facet_wrap(~DoseNum)
# ggplot(CF,aes(y=COMPtot,x=FACIL,col=(time),group=DoseNum))+geom_path()#+facet_wrap(~DoseNum)
#ggplot(CF[time%in%c(21)],aes(y=(COMPtot/FACIL),x=(DoseNum)))+geom_point()+geom_smooth(se=F,method="gam",formula = y~s(x,k=6))
#ggplot(CF[time%in%c(1,6,11,16,21)],aes(y=(FACIL/COMPtot),x=(DoseNum),col=time,group=time))+geom_point()+
 # geom_smooth(se=F,method="gam",formula = y~s(x,k=8))+
 # theme_classic()

#m1<-mgcv::gam(log(I(FACIL/COMPtot))~te(DoseNum,time,k=c(6,8)),data=CF)
#mgcv::vis.gam(m1,plot.type = "contour")
#newY<-expand.grid(DoseNum=seq(0,2000,length=100),time=seq(1,21,length=100))
#newY$Pred<-predict(m1,newdata=newY)
#ggplot(newY[],aes(fill=exp(Pred),y=DoseNum,x=time,group=time))+#geom_point()+
#  geom_raster()+scale_fill_viridis_c(option=2)+theme_classic()+
 # labs(y="Drug dose",x="Time (days)")+
#  scale_y_continuous(breaks=c(0,200,400,600,1000,2000),labels=c(0,200,400,600,1000,2000))
#ggplot(CF[],aes(y=(FACIL/COMPtot),x=(DoseNum),col=time,group=time))+#geom_point()+
#  geom_line(data=newY,aes(y=Pred,x=(DoseNum),col=time,group=time),size=1.5)+
 # theme_classic()
# 
# m1<-mgcv::gam(log(I(FACIL/COMPtot))~te(log(1+DoseNum),time,k=c(6,8)),data=CF)
# #mgcv::vis.gam(m1,plot.type = "contour")
# newY<-expand.grid(DoseNum=  exp( seq(0,7.602,length=1000))-1,time=seq(1,21,length=100))
# newY$Pred<-predict(m1,newdata=newY)
# ggplot(newY[],aes(fill=exp(Pred),y=log(1+DoseNum),x=time,group=time))+#geom_point()+
#   geom_raster(interpolate =T)+scale_fill_viridis_c(option=2)+theme_classic()+
#   labs(y="Drug dose",x="Time (days)")+
#   scale_y_continuous(breaks=log(1+c(0,200,400,600,1000,2000)),labels=c(0,200,400,600,1000,2000))

m1<-mgcv::gam(log(I(FACIL/COMPtot))~te(sqrt(DoseNum),time,k=c(6,8)),data=CF)
#mgcv::vis.gam(m1,plot.type = "contour")
newY<-expand.grid(DoseNum=  ( seq(0,44.73,length=1000))^2,time=seq(1,21,length=100))
newY$Pred<-predict(m1,newdata=newY)

p3<-ggplot(newY[],aes(fill=exp(Pred),y=sqrt(DoseNum),x=time,group=time))+#geom_point()+
  geom_raster(interpolate =T)+scale_fill_viridis_c(name="Facilitation \n relative to \n Competition",limits=c(0,25),option=2)+theme_classic()+
  labs(y="Drug dose",x="Time (days)")+
  scale_y_continuous(breaks=sqrt(c(0,200,400,600,1000,2000)),labels=c(0,200,400,600,1000,2000))+
  theme(aspect.ratio=1)

ggsave(p3,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Facilitation vs competition.png")

p3leg <- ggpubr::get_legend(p3+theme(legend.text = element_blank(),legend.title = element_blank()))
p3LegBlank <- ggpubr::as_ggplot(p3leg)
ggsave(p3LegBlank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Facilitation vs competition Legend.png",height=2,width = 2)

p3Blank <- ggplot(newY[],aes(fill=exp(Pred),y=sqrt(DoseNum),x=time,group=time))+#geom_point()+
  geom_raster(interpolate =T)+scale_fill_viridis_c(name="Facilitation \n relative to \n Competition",limits=c(0,25),option=2)+theme_classic()+
  labs(y="Drug dose",x="Time (days)")+
  scale_y_continuous(breaks=sqrt(c(0,200,400,600,1000,2000)),labels=c(0,200,400,600,1000,2000))+
  theme(aspect.ratio=1)+theme(axis.title = element_blank(),axis.text =  element_blank(),strip.background = element_blank(),
      strip.text.x = element_blank(),legend.position="none")

ggsave(p3Blank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Facilitation vs competition BLANK.png")


m1<-mgcv::gam(log(I(FACIL/COMPtot))~te(sqrt(1+TotN),sqrt(1+TotA),k=c(8,8)),gamma=8,data=CFmono_co)
#mgcv::vis.gam(m1,plot.type = "contour")
newY<-expand.grid(TotN=  ( seq(0,sqrt(max(CFmono_co$TotN)),length=500))^2,
                  TotA=  ( seq(0,sqrt(max(CFmono_co$TotA)),length=500))^2)
newY$Pred<-predict(m1,newdata=newY)
p2<-ggplot(newY[],aes(fill=exp(Pred),y=sqrt(TotN),x=sqrt(TotA)))+#geom_point()+
  geom_raster(interpolate =T)+scale_fill_viridis_c(limits=c(0,25),name="Facilitation \n relative to \n Competition",option=2)+theme_classic()+
  labs(y="Sensitive cell abundance",x="Resistant cell abundance")+
  scale_x_continuous(breaks=sqrt((seq(0,500,by=50))^2),labels=(seq(0,500,by=50))^2)+
  scale_y_continuous(breaks=sqrt((seq(0,500,by=50))^2),labels=(seq(0,500,by=50))^2)+
  theme(aspect.ratio=1)
  
p1
p2


#####
# 
# 
# 
# ggplot(rbind(CF,CF2)[],aes(col=log((FACIL/COMPtot)),y=logit(TotN/(TotN+TotA)),x=sqrt(DoseNum)))+geom_point()+scale_color_viridis_c(option=2)+theme_classic()
# ggplot(rbind(CF,CF2)[],aes(col=log((FACIL/COMPtot)),y=logit(TotN/(TotN+TotA)),x=sqrt(TotN+TotA)))+geom_point()+scale_color_viridis_c(option=2)+theme_classic()
# 
# 
# ggplot(rbind(CF,CF2)[],aes(col=log((FACIL/COMPtot)),y=sqrt(TotN),x=sqrt(TotA)))+geom_point()+scale_color_viridis_c(option=2)+theme_classic()
# 
# 
# ggplot(CF[],aes(y=log((FACIL)),x=sqrt(time), col=sqrt(DoseNum)))+geom_point()+scale_color_viridis_c(option=2)+theme_classic()
# 
# 
#   #geom_point()+
#   geom_raster(interpolate =T)+scale_fill_viridis_c(option=2)+theme_classic()+
#   labs(y="Drug dose",x="Time (days)")+
#   scale_y_continuous(breaks=sqrt(c(0,200,400,600,1000,2000)),labels=c(0,200,400,600,1000,2000))
# 
# CF
# 
# ggplot(CF[time%in%c(1,5,10,15,21)],aes(y=log(COMPtot/FACIL),x=log(1+DoseNum),group=DoseNum))+geom_point()+facet_wrap(~time)
# ggplot(CF[],aes(y=log(COMPtot/FACIL),x=time, col=log(1+DoseNum),group=DoseNum))+geom_path()
