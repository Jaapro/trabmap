trabmap_install_dependencies <- function(){

  if(length(c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr", "tidyverse")[!c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr", "tidyverse") %in% installed.packages()[,"Package"]]) > 0){
    install.packages(c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr", "tidyverse")[!c("EBImage", "tiff", "oro.nifti", "foreach", "mmand", "doParallel", "BiocManager", "iterators", "parallel", "doSNOW", "progress", "plyr", "tidyverse") %in% installed.packages()[,"Package"]])
  }

  remotes::install_github("AlessioVeneziano/IndianaBones")
  BiocManager::install("EBImage")
}
