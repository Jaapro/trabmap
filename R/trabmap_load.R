trabmap_load <- function(){

  if(length(c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr")[!c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr") %in% installed.packages()[,"Package"]]) > 0){
    install.packages(c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr")[!c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr") %in% installed.packages()[,"Package"]])
  }


  library(EBImage)
  library(tiff)
  library(oro.nifti)
  library(foreach)
  library(mmand)
  library(doParallel)
  library(BiocManager)
  library(iterators)
  library(parallel)
  library(doSNOW)
  library(progress)
  library(plyr)
}
