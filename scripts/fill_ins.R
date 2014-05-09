#!/usr/bin/env Rscript
paste("Processing ",commandArgs(TRUE))
r1p<-read.table(commandArgs(TRUE)[1])
r2p<-read.table(commandArgs(TRUE)[2])

logtrans <- function(x){
	y<-log(x)
	y[is.infinite(y)] <- 0
	y
}
logratios<-log(r1p$V2[r1p$V2!=0&r2p$V2!=0])/log(r2p$V2[r1p$V2!=0&r2p$V2!=0])

library("vioplot")
pdf("output/fill_in_read_counts.pdf",width=6.5,height=3.5)
par(cex=0.8)
vioplot(logtrans(r1p$V2/2),logtrans(r2p$V2/2),logratios,col="gray",names=c("log LEnd+Int reads","log REnd+Int reads", "ratio"))
title("Number of fill-in reads per full length template",ylab="log units")
dev.off()
