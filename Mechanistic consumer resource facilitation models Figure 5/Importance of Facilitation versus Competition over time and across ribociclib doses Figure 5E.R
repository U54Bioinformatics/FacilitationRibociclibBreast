### continue on from Facilitation modification simulation script: Facilitation modification fulvestrant predictions simulations Figure 5D.R

# Decompose division rate into positive effect of facilitation and negative effect of competition
facil_comp_isolation <- function(t,Compos ,Dose , r, a, pars1, Y) {
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
    data.table(time= t,
               Compostition= Compos,
               DoseNum= Dose,rho= r,alpha= a,
               COMPA= COMPA, COMPN= COMPN, COMPtot= COMPtot, FACIL= FACIL, CF= CF, r= r_N, TotN= TotN, TotA= TotA)
  } )
}
#facil_comp_isolation(t= 1, Compos= "polyculture" , Dose= 0 , r= 1, a= 1, pars1= pars1, Y= res)

key1 <- expand.grid(time= 1:21, DoseNum= res$DoseNum%>%unique)
CF <- rbindlist( lapply(1:nrow(key1) , function(i){
  tm <- key1[i,]$time
  ds <- key1[i,]$DoseNum
  y <- facil_comp_isolation(t=tm, Compos="polyculture" , Dose= ds , r= 1, a= 1, pars1= pars1, Y= res)
  y$whichSpp <- c("both")
  return(y)
}))

m1 <- mgcv::gam(log(I(FACIL/COMPtot)) ~ te(sqrt(DoseNum), time, k= c(6,8)), data= CF)
newY <- expand.grid(DoseNum=  ( seq(0,44.73,length=1000))^2, time= seq(1,21,length= 100))
newY$Pred <- predict(m1, newdata=newY)
newY
Fig5EModelPredictionData <- newY
#save( Fig5EModelPredictionData, file="/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/Fig5Edata.RData")

p3 <- ggplot(newY[],aes(fill=exp(Pred),y=sqrt(DoseNum),x=time,group=time))+#geom_point()+
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