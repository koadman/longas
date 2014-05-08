#!/usr/bin/env Rscript
paste("Processing ",commandArgs(TRUE))
abc<-read.table(commandArgs(TRUE)[1],header=T,sep="\t")

library("vioplot")
pdf("output/recombinants.pdf",width=6.5,height=3.5)
par(cex=0.8)
vioplot(abc$Count[abc$Recombinant1==0&abc$Recombinant2==0],abc$Count[abc$Recombinant1==1|abc$Recombinant2==1],abc$Max.with.r1[abc$Recombinant1==1],abc$Max.with_r2[abc$Recombinant2==1],names=c("Nonrecombinant","Recombinant","Left progenitors","Right progenitors"),col="gray")
dev.off()
