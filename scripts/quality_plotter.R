#!/usr/bin/env Rscript

library(plotrix)
mylegend <- function(x, y, xlen, ylen, main, tiks,colors){
    rect(x - xlen/10, y-ylen/10,x+xlen + xlen/10,y+ylen-ylen/20, col=rgb(0.95,0.95,0.95,0.8))
    text(x, y+2*ylen/3, main, adj=c(0,0), cex=0.8)
    color.legend(x, y, x+xlen, y+ylen/4, legend=tiks,rect.col=colors, cex=0.7)
}


abc<-read.table("data/full_length/quals.tsv",row.names=1,sep="\t",na.strings="",comment.char="",quote="",header=F)
abc<-as.matrix(abc)
# remove columns which have < 75% of seqs represented
for(i in 1:ncol(abc)){
	if(sum(is.na(abc[,i])) > nrow(abc) / 4){
		abc[,i] <- NA
	}
}

steps <- 100
medd<-apply(abc,2,quantile,probs=seq(0,1,1/steps),na.rm=T)
inds <- which(!is.na(medd[1,]) & !is.na(medd[steps,]))
inds <- inds[inds>40 & inds<1410]


pdf("output/qplot.pdf",width=9,height=4)
plot( x=c(1,1430), y=range(as.numeric(medd),na.rm=T), pch=NA,ylim=c(0,100),main="Base call qualities in 16S sequence",xlab="Position in Infernal bacterial 16S model",ylab="PHRED score" )
colors <- vector()
for(i in 1:(steps/2)){
	cval = 0.9-0.9*(i/(steps/2)) # seems to produce a nice range of constrast for the eye
	colors[i] <- rgb(cval,cval,cval)
	polygon( x=c(inds,inds[length(inds):1]), y=c(medd[i,inds],medd[(steps+2)-i,inds[length(inds):1]]), col=colors[i], border=NA )
}
lines( x=inds, y=medd[(steps/2),inds], lwd=2 )
mylegend(1200,10,200,20,"Quantiles",seq(0,100,by=20),c(colors,rev(colors)))
dev.off()

#
# now plot coverage
#

covcov<-read.table("data/full_length/all_clean_try5.coverages.txt",row.names=1,sep="\t",na.strings="",comment.char="",quote="",header=F)
covcov<-as.matrix(covcov)
# remove columns which have < 75% of seqs represented
for(i in 1:ncol(covcov)){
	if(sum(is.na(covcov[,i])) > nrow(covcov) / 4){
		covcov[,i] <- NA
	}
}

covmed<-apply(covcov,2,quantile,probs=seq(0,1,1/steps),na.rm=T)
covmed<-log(covmed)
pdf("output/covplot.pdf",width=9,height=4)
plot( x=c(1,1430), y=range(as.numeric(covmed),na.rm=T), pch=NA,ylim=c(0,log(6000)),main="Depth of coverage in 16S microassemblies",xlab="Position in Infernal bacterial 16S model",ylab="log Number of reads covering site" )
for(i in 1:(steps/2)){
	indsi <- which(is.finite(covmed[i,]))
	indsi <- indsi[indsi>40 & indsi<1410]
	indsii <- which(is.finite(covmed[steps+2-i,]))
	indsii <- indsii[indsii>40 & indsii<1410]

#	cval = 0.9-(i/(steps/1.3)) # seems to produce a nice range of constrast for the eye
#	colors[i] <- rgb(cval,cval,cval)
	polygon( x=c(rev(indsii),indsi), y=c(covmed[(steps+2)-i,rev(indsii)],covmed[i,indsi]), col=colors[i], border=NA )
}
lines( x=indsi, y=covmed[(steps/2),indsi], lwd=2 )
mylegend(1200,6,200,2.5,"Quantiles",seq(0,100,by=20),c(colors,rev(colors)))

dev.off()

