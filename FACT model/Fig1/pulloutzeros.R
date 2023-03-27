zips <- NULL
zips.rK <- NULL
for (nm in camafiles) {
  tmp <- get(nm)
  tmp <- subset(tmp,level2==0 & level1 <= 200 & env != "M")
  zips <- rbind(zips,tmp)

  tmp.rK <- get(paste0(nm,".rK"))
  tmp.rK <- subset(tmp.rK,level2==0 & level1 <= 200 & env != "M")
  zips.rK <- rbind(zips.rK,tmp.rK)
}

write.csv(zips,"camazips.csv")
write.csv(zips.rK,"camazips.rK.csv")

