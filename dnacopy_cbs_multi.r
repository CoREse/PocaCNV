.libPaths( c( .libPaths(), "R/library") )
if (!requireNamespace("DNAcopy", quietly = TRUE))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org', lib='R/library')
BiocManager::install("DNAcopy",update=FALSE, lib="R/library")
}
library(DNAcopy)

#rddata<-data.frame(mrd=c(1,1,1,1,1,1,5,5,5,5,5,5,5,5,5))
#rddata2<-data.frame(mrd=c(2,2,2,2,2,2,5,5,5,5,5,5,5,5,5))
#rddata$chr <- 1
#rddata$loc <- 1:nrow(rddata)
#rddata2$chr <- 1
#rddata2$loc <- 1:nrow(rddata)

CBSSegment<-function(x)
{
    .libPaths( c( .libPaths(), "R/library") )
    library(DNAcopy)
    x$chr <- 1
    x$loc <- 1:nrow(x)
    set.seed(0)#to avoid randomization
    segment(smooth.CNA(CNA(x[1],x$chr,x$loc,presorted=TRUE)),verbose=0)
}

#rddatalist<-list(rddata,rddata2)
#rddatalist<-vector("list",2)
#rddatalist<-list()
#c(rddatalist,rddata)
#c(rddatalist,rddata2)
#rddatalist[[1]]<-rddata
#rddatalist[[2]]<-rddata2
#dim(rddatalist)

#rddatalist

#ssegdata<-CBSSegment(rddata)
#ssegdata
segAll<-function(rddatalist)
{
    library(parallel)
    cl=makeCluster(ThreadN)

    segdata <- parLapply(cl,rddatalist,CBSSegment)
    #segdata <- lapply(rddatalist,CBSSegment)

    #segendmean <- data.frame(segdata)[c("loc.end", "seg.mean")]

    set.seed(0)#to avoid randomization
    return(segdata)
}

segAll(rddatalist)