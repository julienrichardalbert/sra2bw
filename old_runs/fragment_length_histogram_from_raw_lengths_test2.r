#!/usr/bin/env Rscript

#args, give the following inputs in order:
# fragsizes (in histogram - from lab script)
# output directory
# output label
# plot label

args = commandArgs(trailingOnly=TRUE)
print(args)

dat <- read.table(args[1],sep="\t",header=F)
colnames(dat) <- c("Size","Count")

raw_num <- rep(dat$Size,dat$Count)


xlimits <- c(min(50,quantile(raw_num,0.02)),max(300,quantile(raw_num,0.98)))
ymax <- max(dat$Count)/1000

pdf(paste0(args[2],"/",args[3],"_fragsize_histogram.pdf"))

    plot(dat$Size,dat$Count/1000,xlim=xlimits,ylim=c(0,ymax),type="l",
         las=1,ylab="Counts (x10^3)",xlab="Fragment size (bp)",bty="n",main=args[4])
    abline(v=c(120,150,500),col=rgb(0,0,0,0.5),lty=2)
    text(xlimits[2]*0.8,ymax*0.9,paste0("Mean = ",round(mean(raw_num),1)))
    
dev.off()

write.table(round(mean(raw_num),1),paste0(args[3],"mean.frag.size.txt"),quote=F,col.names=F,row.names=F)
