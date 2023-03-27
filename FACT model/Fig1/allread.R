source("readc.R")  #CAMA-1 files
camafiles <- c("Estradiol1","Fulvestrant3","Raloxifene3","FOHTamoxifen","Letrozole4")
usefiles <- camafiles

# Create files to save the results
anms <- c("level1","level2","type","env")
for (nm in usefiles) {
  tmp <- get(nm)
  tmp$ilvl1 <- (rank(tmp$level1,ties.method="min")-1)/(nrow(tmp)/length(unique(tmp$level1)))+1
  tmp$ilvl2 <- (rank(tmp$level2,ties.method="min")-1)/(nrow(tmp)/length(unique(tmp$level2)))+1
# Reorder columns
  tmp <- tmp[,c("file","day","level1","level2","type","env","ilvl1","ilvl2","rep","cells")]
  assign(nm,tmp)
# Set up to store parameters: Add columns for r, K, alpha
  tmp <- tmp[duplicated(tmp[,anms])==FALSE,]
  tmp <- tmp[,c("file",anms,"ilvl1","ilvl2")]
# Now add on the data
  tmp[,c("r","K","alpha","alphaval")] <- 0
# Rename
  assign(paste0(nm,".rK"),tmp)
}

rm(anms)
