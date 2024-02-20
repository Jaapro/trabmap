#put nifti files in a folder, read all the nifti files, merge them, and export as merged nifti

trabmap_combine_nifti_b <- function(folder, writefolder, corenumber=1, voxelsize = c(0.3)){

  folderlist <- trabmap_folder_list(folder)

  l <- length(folderlist)


  if(l == 2) {
    ra <- readNIfTI(paste(folder, folderlist[1], sep=""))
    dimsra <- dim(ra)
    rb <- readNIfTI(paste(folder, folderlist[2], sep=""))
    dimsrb <- dim(rb)
    result = array(c(ra,rb),dim = c(dimsra[1],dimsra[2],dimsra[3] + dimsrb[3]))
    rm(ra)
    rm(rb)

  }

  if(l == 3) {
    ra <- readNIfTI(paste(folder, folderlist[1], sep=""))
    dimsra <- dim(ra)
    rb <- readNIfTI(paste(folder, folderlist[2], sep=""))
    dimsrb <- dim(rb)
    rc <- readNIfTI(paste(folder, folderlist[3], sep=""))
    dimsrc <- dim(rc)
    result = array(c(ra,rb,rc),dim = c(dimsra[1],dimsra[2],dimsra[3] + dimsrb[3] + dimsrc[3]))
    rm(ra)
    rm(rb)
    rm(rc)

  }

  if(l == 4) {
    ra <- readNIfTI(paste(folder, folderlist[1], sep=""))
    dimsra <- dim(ra)
    rb <- readNIfTI(paste(folder, folderlist[2], sep=""))
    dimsrb <- dim(rb)
    rc <- readNIfTI(paste(folder, folderlist[3], sep=""))
    dimsrc <- dim(rc)
    rc <- readNIfTI(paste(folder, folderlist[4], sep=""))
    dimsrd <- dim(rd)
    result = array(c(ra,rb,rc, rd),dim = c(dimsra[1],dimsra[2], dimsra[3] + dimsrb[3] + dimsrc[3] + dimsrd[3]))
    rm(ra)
    rm(rb)
    rm(rc)
    rm(rd)

  }

  closeAllConnections()
  gc()

  message(paste("saving merged - ", folderlist[i], sep=""))
  trabnifti<-as.nifti(result)
  trabnifti@pixdim[2:4]<-c(voxelsize,voxelsize,voxelsize)
  writeNIfTI(nim=trabnifti,filename=paste(writefolder,"//",folderlist[1],"_merged", sep=""),gzipped=F)
  rm(result)

  gc()
}


