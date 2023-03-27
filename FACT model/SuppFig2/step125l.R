# Pull out the files with temporary names for looping
nm <- "LY2"
tfile <- get(nm)
tfile.rK <- get(paste0(nm,".rK"))
tmp.rK <- tfile.rK
monorows <- (1:nrow(tmp.rK))[substr(tmp.rK$env,1,1) != "M"]
corows <- (1:nrow(tmp.rK))[substr(tmp.rK$env,1,1) == "M"]
# Just loop over monocultures
for (irow in monorows) {
# Extract the data
  tmp0 <- subset(tfile,type ==tmp.rK$type[irow] & env==tmp.rK$env[irow] &
                  level1==tmp.rK$level1[irow] & level2==tmp.rK$level2[irow])
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)

# Find the best fitting r and K
  Kstart <- max(tmp$cells,20000) #200000
  rstart <- 0.2
  tmp.fit <- optim(c(Kstart,rstart),lfit)$par
  tmp.rK$K[irow] <- tmp.fit[1]
  tmp.rK$r[irow] <- tmp.fit[2]
  print(tmp.rK[irow,])

# Plot them up if you want to
  if (plotm==1) {
    plot(cells ~ day,tmp0,pch=19,col=rep)
    Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmp.rK$K[irow],tmp.rK$r[irow],days)
    lines(days,Cpred,lwd=2,col="black")
    title(main=paste(nm,tmp.rK$type[irow],tmp.rK$env[irow]))
    readline('hit return for next plot> ')
  }
}
# Overwrite the monoculture rows
# Add r and K to the coculture rows
tfile.rK[monorows,c("r","K")] <- tmp.rK[monorows,c("r","K")]
tfile.rK[corows,c("r","K")] <- tmp.rK[monorows,c("r","K")]
assign(paste0(nm,".rK"),tfile.rK)

