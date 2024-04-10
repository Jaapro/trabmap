## merge custom fill and bone to trinary

trabmap_custom_fill_to_trinary_nifti <- function(input_folder, output_folder, voxelsize  ){
  folderlist <- list.files(input_folder)

  expfolder <- c(paste(output_folder,"//", folderlist[1],"_custom_filled",  sep = ""))
  dir.create(expfolder, showWarnings = TRUE, recursive = FALSE, mode = "0777")



  stack_1 <- readNIfTI(paste(f_in, folderlist[1], sep=""))
  gc()
  dims <- dim(stack_1)
  stack_1 <- as.integer(stack_1)
  gc()
  stack_1 <- array(stack_1, c(dims[1], dims[2], dims[3]))
  gc()

  stack_2 <- readNIfTI(paste(f_in, folderlist[2], sep=""))
  gc()
  dims <- dim(stack_2)
  stack_2 <- as.integer(stack_2)
  gc()
  stack_2 <- array(stack_2, c(dims[1], dims[2], dims[3]))
  gc()

    closeAllConnections()
  gc()


  # find which is fill and which is bone by calculating size. largest size is fill
  size_stack_1 <- as.numeric(sum(stack_1))
  size_stack_2 <- as.numeric(sum(stack_2))

  if(size_stack_1 > size_stack_2) {
    fill <- stack_1
    rm(stack_1)
    bone <- stack_2
    rm(stack_2)

  } else {
    fill <- stack_2
    rm(stack_2)
    bone <- stack_1
    rm(stack_1)
  }

  #ensure both stacks overlap
  bone <- bone * fill
  trinary_mask <- bone + fill

  ###############
  ###############    calculate some stuff
  ###############
  message("calculating whole bone summary variables")
  ##### calculate BVTV relative to the volume of the whole bone (including cortical volume)
  ##### the sum of these three may be slightly different from 1, I assume due to the cleaning stage

  whole_bone_bvtv <- sum(bone) / sum(fill)
  whole_bone_marrow_space <- 1 - whole_bone_bvtv

  voxel_cube <- voxelsize^3


  total_volume_mm3 <- sum(fill) * voxel_cube
  marrow_volume_mm3 <- whole_bone_marrow_space * total_volume_mm3
  bone_volume_mm3 <- whole_bone_bvtv * total_volume_mm3

  #save calculations
  df <- as.data.frame(folderlist[1])
  df$voxelsize <- voxelsize
  df$voi_diameter_mm <- "N/A"
  df$voi_interval_mm <- "N/A"
  df$strel_diameter <- "N/A"
  df$iters_mask <- "N/A"
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


  write.csv(df, file = paste(expfolder,"//",folderlist[1],"_results.csv", sep=""),row.names = FALSE)

  ## save output
  rm(stack_1)
  rm(stack_2)
  gc()

  dir.create(paste(expfolder,"//tiffstacks", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

  png(file= paste(expfolder,"//tiffstacks//" ,folderlist[1],"_trinary_xy.png", sep=""), width=350, height=350)
  image(trinary_mask[,,0.5*dims[3]], main = paste(folderlist[1]," trinary XY"))
  dev.off()

  png(file= paste(expfolder,"//tiffstacks//" ,folderlist[1],"_trinary_xz.png", sep=""), width=350, height=350)
  image(trinary_mask[,0.5*dims[2],], main = paste(folderlist[1]," trinary XZ"))
  dev.off()

  png(file= paste(expfolder,"//tiffstacks//" ,folderlist[1],"_trinary_yz.png", sep=""), width=350, height=350)
  image(trinary_mask[0.5*dims[1],,], main = paste(folderlist[1]," trinary YZ"))
  dev.off()
  gc()



  print(paste("Exporting as nifti - trinary mask of : ", folderlist[1]))
  trabnifti<-as.nifti(trinary_mask)
  trabnifti@pixdim[2:4]<-c(voxelsize,voxelsize,voxelsize)
  writeNIfTI(nim=trabnifti,filename=paste(expfolder,"//tiffstacks//",folderlist[1],"_trinary", sep=""),gzipped=F)
  gc()
  rm(trinary_mask)
  rm(trabnifti)
  closeAllConnections()



}

