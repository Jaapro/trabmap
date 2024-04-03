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
    print(paste("calculating dims of: ", folderlist[i]))
    dims <- dim(vol)
    print(paste("vol as int: ", folderlist[i]))
    vol <- as.integer(vol)
    print(paste("vol as int array: ", folderlist[i]))
    vol <- array(vol, c(dims[1], dims[2], dims[3]))
    #vol<-readTiffStack(path=paste(folder, folderlist[i], sep=""), cores=1)
    #memory.limit()
    gc() #clears memory





    #pad the dataset to avoid edge effects later
    #padding the original image stack with zeros so that the edges can be measured with VOIs
    big_array <- array(0, c(dims[1]+2*pad, dims[2] + 2*pad , dims[3])) #creates array filled with 0 that is bigger than the original image based on the voi size
    ba_dims <- dim(big_array)
    big_array <- as.integer(big_array)
    big_array <- array(big_array, c(ba_dims[1], ba_dims[2], ba_dims[3]))
    #create padded image by summing the big array with the smaller vol array
    big_array[(pad + 1):(dims[1] + pad) ,(pad + 1):(dims[2] + pad) , ] <- vol[,,] + big_array[(pad +1 ):(dims[1] + pad) ,(pad + 1):(dims[2] + pad) ,]

    vol<- big_array

    rm(big_array)
    gc()

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
      #mask<-EBImage::closing(mask,d_e_strel)
      mask<-EBImage::fillHull(mask)
      for(j in 1:d_e_steps){mask<-EBImage::erode(mask,d_e_strel)}
    }


    trinary_mask <- vol[,,] + mask #calculate with uncropped data

    if(pad > 0){
      trinary_mask <- trinary_mask[(pad+1):(dims[1]-pad),
                                   (pad+1):(dims[2]-pad),]

      mask <- mask[(pad+1):(dims[1]-pad),
                   (pad+1):(dims[2]-pad),]
    }

    dims <- dim(mask)

    ###############
    ###############    calculate some stuff
    ###############
    message("calculating whole bone summary variables")
    ##### calculate BVTV relative to the volume of the whole bone (including cortical volume)
    ##### the sum of these three may be slightly different from 1, I assume due to the cleaning stage

    whole_bone_bvtv <- sum(vol) / sum(mask)
    whole_bone_marrow_space <- 1 - whole_bone_bvtv

    voxel_cube <- voxelsize^3


    total_volume_mm3 <- sum(mask) * voxel_cube
    marrow_volume_mm3 <- whole_bone_marrow_space * total_volume_mm3
    bone_volume_mm3 <- whole_bone_bvtv * total_volume_mm3

    #save calculations
    df <- as.data.frame(folderlist[i])
    df$voxelsize <- voxelsize
    df$voi_diameter_mm <- "N/A"
    df$voi_interval_mm <- "N/A"
    df$strel_diameter <- strel
    df$iters_mask <- steps
    df$iters_fill <- "N/A"
    df$iters_compact <- "N/A"
    df$clean_voids <- "N/A"
    df$trabecular_bvtv <- "N/A"
    df$whole_bone_bvtv <- whole_bone_bvtv
    df$whole_bone_trab_bvtv <- "N/A"
    df$whole_bone_cort_bvtv <- "N/A"
    df$whole_bone_marrow_bvtv <- whole_bone_marrow_space
    df$total_volume_mm3 <- floor(total_volume_mm3)
    df$trabecular_volume_mm3 <- "N/A"
    df$cortical_volume_mm3 <- "N/A"
    df$marrow_trabecular_volume_mm3 <- "N/A"
    df$marrow_volume_mm3 <- floor(marrow_volume_mm3)
    df$whole_bone_volume_mm3 <- floor(bone_volume_mm3)


    write.csv(df, file = paste(expfolder,"//",folderlist[i],"_results.csv", sep=""),row.names = FALSE)


    rm(vol)

    dir.create(paste(expfolder,"//tiffstacks", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

    png(file= paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill_xy.png", sep=""), width=350, height=350)
    image(mask[,,0.5*dims[3]], main = paste(folderlist[i]," Fill XY"))
    dev.off()

    png(file= paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill_xz.png", sep=""), width=350, height=350)
    image(mask[,0.5*dims[2],], main = paste(folderlist[i]," Fill XZ"))
    dev.off()

    png(file= paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill_yz.png", sep=""), width=350, height=350)
    image(mask[0.5*dims[1],,], main = paste(folderlist[i]," Fill YZ"))
    dev.off()

    png(file= paste(expfolder,"//tiffstacks//" ,folderlist[i],"_trinary_xy.png", sep=""), width=350, height=350)
    image(trinary_mask[,,0.5*dims[3]], main = paste(folderlist[i]," trinary XY"))
    dev.off()

    png(file= paste(expfolder,"//tiffstacks//" ,folderlist[i],"_trinary_xz.png", sep=""), width=350, height=350)
    image(trinary_mask[,0.5*dims[2],], main = paste(folderlist[i]," trinary XZ"))
    dev.off()

    png(file= paste(expfolder,"//tiffstacks//" ,folderlist[i],"_trinary_yz.png", sep=""), width=350, height=350)
    image(trinary_mask[0.5*dims[1],,], main = paste(folderlist[i]," trinary YZ"))
    dev.off()
    gc()


    print(paste("Exporting fill of : ", folderlist[i]))
    dir.create(paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    writeTiffStack(Stack=mask,path=paste(expfolder,"//tiffstacks//" ,folderlist[i],"_fill", sep=""),leadname=paste(folderlist[i],"_fill", sep=""))
    rm(mask)
    gc()


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
