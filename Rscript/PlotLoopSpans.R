### ChIA-Drop stat plot Fragment to Fragment distance
suppressMessages({
  if (!require("ggplot2")) install.packages('ggplot2', repos = "https://cloud.r-project.org/")
  library("ggplot2")})

args = commandArgs(trailingOnly = TRUE)

FIN <- read.table(args[1])

names(FIN) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "len", "fragcount")

MX=round(density(FIN$len)$x[which.max(density(FIN$len)$y)])
MXPOS = paste0("max(density)$x = ", MX)
f2si <- c("500b","10k", "1M", "10M", "100M")
MFN=15

TITLE <- "LoopSpans.dist"
PLOTNAME <- paste0(TITLE, ".png")


png(PLOTNAME, 1200, 1200)
ggplot(FIN,aes(len, color=fragcount,group=fragcount))+
  geom_density(aes(y=..scaled..),size=3)+
  theme_bw(base_size=40)+
  theme(axis.text = element_text(size = 55))+
  theme(legend.title = element_text(size = 40))+
  theme(legend.key.size = unit(0.5, "inches"), legend.text = element_text(size = 30))+
  theme(strip.text = element_text(size = 40))+
  theme(plot.title = element_text(hjust = 0.5, size = 50))+
  theme(axis.title = element_text(hjust = 0.5, size = 75))+
  xlab("log10 Loop spans (bp)")+
  ggtitle(paste0(TITLE))+
  ylab("scaled density")+
  #  scale_x_continuous(limits=c(0,10000000),breaks=c(100,10000,100000,10000000,1000000,MXPOS))+  
  #scale_x_log10(labels=f2si,breaks=c(0,1,100,10000,100000,10000000,1000000,100000000,1000000000))+
  scale_x_log10(breaks=c(500,10^c(4, 6, 7, 8)),labels=f2si,limits=c(500,1000000000))+
  #  scale_colour_gradientn(colours = rev(rainbow(3)),name="FR#")+
  scale_colour_gradientn(colours = rev(rainbow(3)),name="FN",breaks=c(1,round(MFN/2),MFN),labels=c(1,round(MFN/2),MFN),limits=c(1,MFN),na.value = "transparent")+
  theme(panel.grid = element_blank(),
        legend.position="bottom",
        plot.title = element_text(size = 15, face = "bold")
  )

dev.off()
