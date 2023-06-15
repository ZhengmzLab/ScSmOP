# scatter plot for Loop span and pet count
library(dplyr)
library(ggplot2)

RNAP2W5FUT00 <- read.table("E:/P53AnalysisByJK/ChIA-PET/Loop/RNAP2W5FUT00_hg38.e500.cluster")
RNAP2W5FUT09 <- read.table("E:/P53AnalysisByJK/ChIA-PET/Loop/RNAP2W5FUT09_hg38.e500.cluster")

 
All.Loop <- rbind(RNAP2W5FUT00, RNAP2W5FUT09)
All.Loop.Intra <- All.Loop %>% filter(V1 == V4)


All.Loop.Intra[, "LoopSpan"] <- floor((All.Loop.Intra[, 5] + All.Loop.Intra[, 6])/2 - (All.Loop.Intra[, 3] + All.Loop.Intra[, 2])/2)
names(All.Loop.Intra)[7] <- "PETCount"

All.Loop.Intra[, "Log10LoopSpan"] <- log10(All.Loop.Intra$LoopSpan)
All.Loop.Intra[, "Log10PETCount"] <- log10(All.Loop.Intra$PETCount)

p.data <- All.Loop.Intra

p <- ggplot(p.data, aes(x = LoopSpan, y = PETCount)) + geom_point(size = 3, color = "blue", shape = 16)
# p <- ggplot(p.data, aes(x = Rank, y = Log2WTFC, color = KOColor)) + geom_point(size = 3)
# p <- p + scale_color_gradient(low = "white", high = "black")
# p <- p + scale_color_manual(breaks = c(0:10), values = c("white", "black"))
p <- p + scale_x_log10(breaks = c(500,10^c(4, 6, 7, 8)),labels=c("500b","10k", "1M", "10M", "100M"))
p <- p + scale_y_log10(breaks = c(10, 100))
# p <- p + scale_x_continuous(breaks = c(3:9), labels = c("1k","10k", "100k", "1M", "10M", "100M", "1B"))
# p <- p + scale_y_continuous(breaks = c(1,2), labels = c("10", "100"))
p <- p + xlab("Loop span (bp)")
p <- p + ylab("PET Count")
p <- p + ggtitle("RNAP2WT00_RNAP2WT09_LoopSpan-PETCount")
p <- p + theme(axis.text = element_text(size = 55))
p <- p + theme(legend.title = element_text(size = 40))
p <- p + theme(legend.key.size = unit(0.5, "inches"), legend.text = element_text(size = 30))
p <- p + theme(strip.text = element_text(size = 40))
p <- p + theme(plot.title = element_text(hjust = 0.5, size = 50))
p <- p + theme(axis.title = element_text(hjust = 0.5, size = 75))
p <- p + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", linetype = "solid", size = 2, fill = NA))

png("E:/P53AnalysisByJK/ChIA-PET/RNAP2WT00_RNAP2WT09_LoopSpan-PETCount.png", 1400, 1200)
print(p)
dev.off()


FIN <- All.Loop.Intra[, c(1,8,7)]

names(FIN) <- c("name", "len", "fragcount")

MX=round(density(FIN$len)$x[which.max(density(FIN$len)$y)])
MXPOS = paste0("max(density)$x = ", MX)
f2si <- c("500b","10k", "1M", "10M", "100M")
MFN=15

TITLE <- "E:/P53AnalysisByJK/ChIA-PET/RNAP2WT00_RNAP2WT09_LoopSpan-PETCount_Density"
PLOTNAME <- paste0(TITLE, ".png")


png(PLOTNAME, 1400, 1200)
p <- ggplot(FIN,aes(len)) + geom_density(aes(y=..scaled..),size=3)
p <- p + xlab("Loop span (bp)")
p <- p + ggtitle(paste0(TITLE))
p <- p + ylab("scaled density")
p <- p + scale_x_log10(breaks=c(500,10^c(4, 6, 7, 8)),labels=f2si) 
p <- p + scale_colour_gradientn(colours = rev(rainbow(3)),name="FN",breaks=c(1,round(MFN/2),MFN),labels=c(1,round(MFN/2),MFN),limits=c(1,MFN),na.value = "transparent")
p <- p + theme(panel.grid = element_blank(),legend.position="bottom", plot.title = element_text(size = 15, face = "bold"))
p <- p + theme(axis.text = element_text(size = 55))
p <- p + theme(legend.title = element_text(size = 40))
p <- p + theme(legend.key.size = unit(0.5, "inches"), legend.text = element_text(size = 30))
p <- p + theme(strip.text = element_text(size = 40))
p <- p + theme(plot.title = element_text(hjust = 0.5, size = 50))
p <- p + theme(axis.title = element_text(hjust = 0.5, size = 75))
p <- p + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", linetype = "solid", size = 2, fill = NA))
print(p)
dev.off()
