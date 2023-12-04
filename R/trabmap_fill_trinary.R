################
#### Function to fill bone images and produce trinary masks (e.g. bone, marrow, air)
#### Trinary masks serve as input for creating 3D heatmaps of bone volume fraction
#################

trabmap_fill_trinary <- function(folder, writefolder, voxelsizelist = c(), strel = 11, steps = 4, pad=TRUE,df_list = list()) {
  d_e_strel <-  makeBrush(size=strel,shape="disc") # size of filling element in pixels
  d_e_steps <-  steps  # number of initial erosion and dilation steps to close the shell
  folderlist <- list.files(folder)

  #padding to avoid potential edge effects
  if (pad == TRUE){
    pad <- steps*strel
  } else{
    pad <- 0
  }


  for (i in 1:length(folderlist)) {
    voxelsize <- voxelsizelist[i]


    print(paste("now running: ", folderlist[i]))
    vol<-trabmap_readtiffstack(path=paste(folder, folderlist[i], sep=""), cores=1)
    #vol<-readTiffStack(path=paste(folder, folderlist[i], sep=""), cores=1)
    #memory.limit()
    gc() #clears memory


    dims <- dim(vol)


    #pad the dataset to avoid edge effects later
    #padding the original image stack with zeros so that the edges can be measured with VOIs
    big_array <- array(0, c(dims[1]+2*pad, dims[2] + 2*pad ,dims[3] + 2*pad)) #creates array filled with 0 that is bigger than the original image based on the voi size
    #create padded image by summing the big array with the smaller vol array
    big_array[(pad + 1):(dims[1] + pad) ,(pad + 1):(dims[2] + pad) , (pad +1):(dims[3] + pad)] <- vol[,,] + big_array[(pad +1 ):(dims[1] + pad) ,(pad + 1):(dims[2] + pad) ,(pad +1):(dims[3] + pad)]

    vol<- big_array
    rm(big_array)

    dims <- dim(vol)
    n <- 0.5*dims[3]
    #export files

    expfolder <- c(paste(writefolder,"//", folderlist[i],"_IB",  sep = ""))
    dir.create(expfolder, showWarnings = TRUE, recursive = FALSE, mode = "0777")


    message("now performing: dilation, erosion, and fill")
    mask<-vol
    if(length(d_e_steps)==0){
      mask<-EBImage::fillHull(mask)
    } else {
      for(j in 1:d_e_steps){mask<-EBImage::dilate(mask,d_e_strel)}
      mask<-EBImage::closing(mask,d_e_strel)
      for(j in 1:d_e_steps){mask<-EBImage::erode(mask,d_e_strel)}
      mask<-EBImage::fillHull(mask)
    }


    trinary_mask <- vol[,,] + mask #calculate with uncropped data

    if(pad > 0){
      trinary_mask <- trinary_mask[(pad+1):(dims[1]-pad),
                                   (pad+1):(dims[2]-pad),
                                   (pad+1):(dims[3]-pad)]

      mask <- mask[(pad+1):(dims[1]-pad),
                   (pad+1):(dims[2]-pad),
                   (pad+1):(dims[3]-pad)]
    }




    dir.create(paste(expfolder,"//tiffstacks", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

    print(paste("Exporting fill of : ", folderlist[i]))
    dir.create(paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    writeTiffStack(Stack=mask,path=paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill", sep=""),leadname=paste(folderlist[i],"_fill", sep=""))
    rm(mask)


    print(paste("Exporting as nifti - trinary mask of : ", folderlist[i]))
    dir.create(paste(expfolder,"//tiffstacks//" ,folderlist[i],"_trinary_mask", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    trabnifti<-as.nifti(trinary_mask)
    trabnifti@pixdim[2:4]<-c(voxelsize,voxelsize,voxelsize)
    writeNIfTI(nim=trabnifti,filename=paste(expfolder,"//tiffstacks//",folderlist[i],"_trinary_mask//",folderlist[i],"_trinary_mask", sep=""),gzipped=F)
    gc()
    rm(trinary_mask)
    rm(trabnifti)

    closeAllConnections()

  }


}
