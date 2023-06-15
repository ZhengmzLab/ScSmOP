### ChIA-Drop stat plot Fragment Length
suppressMessages({
  if (!require("ggplot2")) install.packages('ggplot2', repos = "https://cloud.r-project.org/")
  library("ggplot2")})


args = commandArgs(trailingOnly = TRUE)

FIN <- read.table(args[1])

names(FIN) <- c("chr", "start", "end", "name", "len")

MX = round(density(FIN$len)$x[which.max(density(FIN$len)$y)])
MXPOS = paste0("max(density)$x = ", MX)


TITLE <- "FragmentLength.size"
PLOTNAME <- paste0(TITLE, ".png")
png(PLOTNAME, 400, 400)


ggplot(FIN,aes(len))+
  geom_density(aes(y=..scaled..),size=1,color="navy")+
  theme_bw(base_size=20)+
  xlab("fragment size (bp)")+
  ggtitle(paste0(TITLE," \n ",MXPOS))+
  ylab("scaled density")+
  #ggtitle(LIB61)+
  scale_x_continuous(limits=c(0,800))+
  
  theme(panel.grid = element_blank(),
        legend.position="none",
        plot.title = element_text(size = 15, face = "bold")
  )
dev.off()
