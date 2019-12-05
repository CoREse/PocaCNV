if (!requireNamespace("DNAcopy", quietly = TRUE))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
BiocManager::install("DNAcopy",update=FALSE)
}
library(DNAcopy)

rddata

rddata$chr <- 1
rddata$loc <- 1:nrow(rddata)

segdata <- segment(CNA(rddata[1],rddata$chr,rddata$loc),verbose=0)

#segendmean <- data.frame(segdata)[c("loc.end", "seg.mean")]

segdata
