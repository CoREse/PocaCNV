.libPaths( c( .libPaths(), "R/library") )
if (!requireNamespace("DNAcopy", quietly = TRUE))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org', lib='R/library')
BiocManager::install("DNAcopy",update=FALSE, lib="R/library")
}
library(DNAcopy)

rddata$chr <- 1
rddata$loc <- 1:nrow(rddata)

segdata <- segment(smooth.CNA(CNA(rddata[1],rddata$chr,rddata$loc,presorted=TRUE)),verbose=0)

#segendmean <- data.frame(segdata)[c("loc.end", "seg.mean")]

segdata
