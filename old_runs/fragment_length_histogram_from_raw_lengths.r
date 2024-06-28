#!/usr/bin/env Rscript

#read in fragment length file name and output prefix
   #args <- c("King2017_Brg1_rep1.fraglengths.txt","King2017_Brg1_rep1")
   args = commandArgs(trailingOnly=TRUE)
   print(args)
   #1) Frag lengths
   #2) output prefix
   
#read in fragment lengths
   dat <- read.table(args[1])$V1
#remove lengths <= 0 or > 1000, typically this is 1/100000
   print(paste0("removing ",sum(!(dat > 0 & dat <= 1000))," fragment lengths <= 0 or > 1000"))
   print(paste0("now calculating histogram on ",sum(dat > 0 & dat <= 1000)," remaining fragment lengths"))
   dat <- dat[dat > 0 & dat <= 1000]

#Calculate histogram and make dataframe with values to plot and output in text file
   dat_hist <- hist(dat,breaks=seq(0.5,1000.5,1),plot=F)
   dat_df <- data.frame(Fragment_size=dat_hist$mids,
                        counts=dat_hist$counts,
                        counts_per_mil=dat_hist$counts / (length(dat)/1000000),
                        density=dat_hist$density)
#Write histogram value to txt
   write.table(dat_df,paste0(args[2],"_fragment_size_histogram.txt"),sep="\t",quote=F,row.name=F)

#calculate peak size
   peakSize <- dat_df$Fragment_size[which.max(dat_df$counts_per_mil)]

#output fragment size summary stats to txt file
   write.table(data.frame(c(names(summary(dat)),"PeakSize"),
                          c(as.numeric(summary(dat)),peakSize)),
               paste0(args[2],"_fragment_size_summary_stats.txt"),sep="\t",quote=F,
               row.names=F,col.names=F)
#make PDF of histogram counts / million fragments, plot x-axis from 0.05 to 99.5 % values, draw a vertical line at the peak frag size   
   xlimits <- quantile(dat,c(0.005,0.995))
   pdf(paste0(args[2],"_fragment_size_histogram_counts_per_million.pdf"),height=4,width=4,useDingbats=F)
   plot(dat_df$Fragment_size,
        dat_df$counts_per_mil,
        type="l",xlim=xlimits,las=1,
        ylab="Counts per million",
        xlab="Fragment length (bp)",bty="n",
        main=paste0(args[2],"(",round(length(dat)/1000000,2),")"),cex=0.9)
   legend("topright",legend=paste0("Peak size = ",peakSize),bty="n")
   abline(v=peakSize,lty=2,col=rgb(0,0,0,0.5))
   dev.off()
   
