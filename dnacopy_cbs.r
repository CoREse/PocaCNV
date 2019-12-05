if (!requireNamespace("DNAcopy", quietly = TRUE))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",quiet=TRUE)
BiocManager::install("DNAcopy")
}
library(DNAcopy)

rddata

rddata$chr <- 1
rddata$loc <- 1:nrow(rddata)

segdata <- segment(CNA(rddata[1],rddata$chr,rddata$loc),verbose=0)

#segendmean <- data.frame(segdata)[c("loc.end", "seg.mean")]

segdata
