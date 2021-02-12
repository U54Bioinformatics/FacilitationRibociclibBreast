rm(list=ls())
require(deSolve) ;require(ggplot2) ; require(tidyr) ; require(dplyr) ; require(data.table) 
require(Matrix); require(fda); library("readr"); require(boot);
require(stringr)
require(loo)
setwd("/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/")
filenms0<-list.files()
filenms<-filenms0[][!grepl(pattern="old",filenms0)&grepl("Ribo",filenms0)]

com<-rbindlist(lapply(1:length(filenms),function(x){
  cat(x/length(filenms));cat("     ");
  load(file=filenms[x])
  log_lik_1=extract_log_lik(test, merge_chains = F)
  r_eff_1=relative_eff(log_lik_1)
  loo_1 <- waic(log_lik_1)
  res<-data.table(as.data.table(loo_1$ estimates,keep.rownames = T),mod=filenms[x])
  return(res)#sum(loo_1$ pointwise[,"waic"])
}))
com[,ucl:=Estimate+1.96*SE]
com[,lcl:=Estimate-1.96*SE]
com[,Model:="old" ]

com[mod=="Ribociclib facilitation life history FULL NQS.Rdata",Model:="Facilitation (LH1)" ]
com[mod=="Ribociclib facilitation life history Simplified NS.Rdata",Model:="Facilitation (LH2)" ]
com[mod=="Ribociclib facilitation happy cell N.Rdata",Model:="Happy cells" ]
com[mod=="Ribociclib facilitation Allee effect N.Rdata",Model:="Allee effect" ]
com[mod=="Ribociclib Plasticity random switching.Rdata",Model:="Plasticity (GC)" ]
com[mod=="Ribociclib Plasticity drug induced switching.Rdata",Model:="Plasticity (DI)" ]
com[mod=="Ribociclib subclonal competition N.Rdata",Model:="Subclonal evolution" ]
com$Model <- factor(com$Model, levels = rev(c("Facilitation (LH1)","Facilitation (LH2)","Happy cells","Allee effect","Plasticity (GC)",
                                          "Plasticity (DI)","Subclonal evolution")))

com[,Modelgroup:="Facilitation"]

com[Model%in%c("Plasticity (DI)","Plasticity (GC)"),Modelgroup:="Plasticity"]
com[Model%in%c("Subclonal evolution"),Modelgroup:="Selection"]

p4<-ggplot(com[rn=="waic"][Model!="old"] ,aes(x=Estimate,y=Model,col=Modelgroup))+geom_point(size=6)+
  geom_errorbarh(aes(xmin=lcl,xmax=ucl),height=0.4,size=1.3)+
  theme_classic(base_size=21)+labs(x="Prediction error (WAIC)")+theme(aspect.ratio=1)+
  scale_color_manual(name="Process",values=RColorBrewer::brewer.pal(3, "Set2"))

ggsave(p4,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Model comparisons.png")

p4leg <- ggpubr::get_legend(p4+theme(legend.text = element_blank(),legend.title = element_blank()))
p4LegBlank <- ggpubr::as_ggplot(p4leg)
ggsave(p4LegBlank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Model comparisons Legend.png",height=2,width = 2)

p4Blank <- ggplot(com[rn=="waic"][Model!="old"] ,aes(x=Estimate,y=Model,col=Modelgroup))+geom_point(size=6)+
  geom_errorbarh(aes(xmin=lcl,xmax=ucl),height=0.4,size=1.3)+
  theme_classic(base_size=21)+labs(x="Prediction error (WAIC)")+theme(aspect.ratio=1)+
  scale_color_manual(name="Process",values=RColorBrewer::brewer.pal(3, "Set2"))+
  theme(aspect.ratio=1)+theme(axis.title = element_blank(),axis.text =  element_blank(),strip.background = element_blank(),
                              strip.text.x = element_blank(),legend.position="none")

ggsave(p4Blank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Model comparisons BLANK.png")



