suppressMessages({
  if (!require("DropletUtils"))
  {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE)
    BiocManager::install("DropletUtils", update = TRUE)
  }
  library("Matrix")
  library("DropletUtils")
})

print("Reading peaks and cells")
barcodes <- read.table("Cell")
feature <- read.table("Peak")
# for(i in seq_len(nrow(feature)))
# {
#   feature[i, 4] <- paste0("peak_", i)
# }
print("Reading matrix")
mat <- readMM("matrix.mtx")


colnames(mat) <- barcodes$V1
rownames(mat) <- feature$V4
print("Filtering cells")
bcrank <- barcodeRanks(mat)
uniq <- !duplicated(bcrank$rank) 

OutputDrop <- emptyDrops(mat)

FinalDrop <- OutputDrop[which(OutputDrop$FDR<=0.001), ]

write.table(as.data.frame(row.names(FinalDrop)), "filtered_cells.tsv", row.names = F, quote = F, col.names = F)

