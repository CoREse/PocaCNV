if (!requireNamespace("DNAcopy", quietly = TRUE))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
BiocManager::install("DNAcopy",update=FALSE)
}

