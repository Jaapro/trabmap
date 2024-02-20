#put nifti files in a folder, read all the nifti files, merge them, and export as merged nifti

trabmap_combine_nifti <- function(folder, writefolder, corenumber=1, voxelsize = c(0.3)){

  folderlist <- trabmap_folder_list(folder)

  stacklist <- list()



  for (i in 1:length(folderlist)) {
    message(paste("running whole bone analysis - ", folderlist[i], sep=""))

    gc()

    stacklist[[i]] <- readNIfTI(paste(folder, folderlist[i], sep=""))

  }

  l <- length(folderlist)


  if(l == 2) {
    ra<-simplify2array(stacklist[[1]])
    dims <- dim(ra)
    rb <- simplify2array(stacklist[[2]])
    result = array(c(ra,rb),dim = c(dims[1],dims[2],dims[3]*l))

  }

  if(l == 3) {
    ra<-simplify2array(stacklist[[1]])
    dims <- dim(ra)
    rb <- simplify2array(stacklist[[2]])
    rc <- simplify2array(stacklist[[2]])
    result = array(c(ra,rb, rc),dim = c(dims[1],dims[2],dims[3]*l))

  }

  if(l == 4) {
    ra<-simplify2array(stacklist[[1]])
    dims <- dim(ra)
    rb <- simplify2array(stacklist[[2]])
    rc <- simplify2array(stacklist[[3]])
    rd <- simplify2array(stacklist[[4]])
    result = array(c(ra, rb, rc, rd),dim = c(dims[1],dims[2],dims[3]*l))

  }

  if(l == 5) {
    ra<-simplify2array(stacklist[[1]])
    dims <- dim(ra)
    rb <- simplify2array(stacklist[[2]])
    rc <- simplify2array(stacklist[[3]])
    rd <- simplify2array(stacklist[[4]])
    re <- simplify2array(stacklist[[5]])
    result = array(c(ra, rb, rc, rd, re),dim = c(dims[1],dims[2],dims[3]*l))

  }
  rm(stacklist)
  rm(ra)
  rm(rb)
  rm(rc)
  rm(rd)
  rm(re)
  closeAllConnections()
  gc()

  message(paste("saving merged - ", folderlist[i], sep=""))
  trabnifti<-as.nifti(result)
  trabnifti@pixdim[2:4]<-c(voxelsize,voxelsize,voxelsize)
  writeNIfTI(nim=trabnifti,filename=paste(writefolder,"//",folderlist[i],"_merged", sep=""),gzipped=F)
  rm(result)

  gc()
}


