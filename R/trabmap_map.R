trabmap_map <- function(folder, writefolder, voxelsizelist = c(), voi_diameter_mm=5, voi_interval_mm=2.5, corenumber=1){

  folderlist <- trabmap_folder_list(folder)

  for (i in 1:length(folderlist)) {

    #
    #
    #
    ##########
    ## create bvtv maps
    #########
    #
    #
    #
    message(paste("running whole bone analysis - ", folderlist[i], sep=""))

    gc()
    #trinary trab mask
    tb_array <- readNIfTI(paste(folder, folderlist[i], sep=""))

    expfolder <- c(paste(writefolder,"//", folderlist[i],"_BVTV",  sep = ""))
    expfolder <- sub("\\.nii", "", expfolder)
    dir.create(expfolder, showWarnings = TRUE, recursive = FALSE, mode = "0777")

    name <- folderlist[i]
    name <- sub("\\.nii", "", name)


    #### make a mesh for interpolation later

    ###############3
    ######## Export ply mesh
    ################

    meshpath <- c(paste(expfolder, "//", name, "_mesh.ply", sep=""))

    message(paste("creating mesh of - ", name, sep=""))
    trabmap_stack_to_ply(image_stack = tb_array,
                         voxelsize = voxelsizelist[i],
                         output_ply_path = meshpath,
                         pad_thickness = 10)
    gc()


    dims <- dim(tb_array)

    tb_array <- as.integer(tb_array)
    gc()

    tb_array <- array(tb_array, c(dims[1], dims[2], dims[3]))


    voxelsize_original <- voxelsizelist[i] # voxelsize of the original image

    voisize <- floor(voi_diameter_mm/voxelsize_original) # VOI diameter in pixels -- MUST BE ODD



    #voisize has to be an odd number to create a circle later, the following makes sure its odd
    if((voisize %% 2)==0) {
      voisize <- voisize +1
    }

    voi_interval <- floor(voi_interval_mm/voxelsize_original) #interval between VOIs in pixels
    voxelsize <- voxelsize_original * voi_interval   #new voxel size of the heatmap imagestack

    #make spherical VOI --- the following loop creates a cubic VOI containing a spherical VOI
    empty_cube <- array(0,c(voisize, voisize, voisize))
    loop_length <- ceiling(voisize/2)
    expand_factor <- 0



    for (l in 1:loop_length){
      if (l != loop_length){
        empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor) , l] <- empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor) , l] + makeBrush(size=l+expand_factor, shape="disc") #start of stack
        empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor)  , voisize - expand_factor] <- empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor) , voisize - expand_factor] + makeBrush(size=l+expand_factor, shape="disc") #end of stack
        expand_factor <- expand_factor + 1
      }
      else {
        empty_cube[ , , loop_length] <- empty_cube[ , , loop_length] + makeBrush(size=voisize, shape="disc") #start of stack
      }
    }


    #padding the original image stack with zeros so that the edges can be measured with VOIs
    dims_slice <- dims
    big_array <- array(0, c(dims_slice[1]+voisize, dims_slice[2] + voisize ,dims_slice[3] + voisize)) #creates array filled with 0 that is bigger than the original image based on the voi size
    ba_dims <- dim(big_array)
    big_array <- as.integer(big_array)
    big_array <- array(big_array, c(ba_dims[1], ba_dims[2], ba_dims[3]))


    #create padded image by summing the big array with the smaller tb_array array
    big_array[(0.5*voisize + 1):(dims_slice[1] + 0.5*voisize) ,(0.5*voisize + 1):(dims_slice[2] + 0.5*voisize) , (0.5*voisize +1):(dims_slice[3] + 0.5*voisize)] <- tb_array[,,] + big_array[(0.5*voisize +1 ):(dims_slice[1] + 0.5*voisize) ,(0.5*voisize + 1):(dims_slice[2] + 0.5*voisize) ,(0.5*voisize +1):(dims_slice[3] + 0.5*voisize)]
    dims_padded <- dim(big_array) #get the length width and depth of the imagestack


    #set the boundaries of the volume within which vois will be taken
    firstpixel_i <- ceiling(0.5*voisize)
    firstpixel_j <- ceiling(0.5*voisize)
    lastpixel_i <- dims_padded[1]-floor(0.5*voisize)
    lastpixel_j <- dims_padded[2]-floor(0.5*voisize)
    firstpixel_z <- ceiling(0.5*voisize)
    lastpixel_z <- as.numeric(dims_padded[3]-floor(0.5*voisize))

    #determine intervals where VOIs will be placed
    loopseq_i <- seq(from = firstpixel_i, to=lastpixel_i, by=voi_interval)
    loopseq_j <- seq(from = firstpixel_j, to=lastpixel_j, by=voi_interval)
    loopseq_z <- seq(from = firstpixel_z, to=lastpixel_z, by=voi_interval)

    #create an empty imagestack that will be filled with new trabecular values
    bvtv_map <- array(0, c(length(loopseq_i),length(loopseq_j),length(loopseq_z)))




    ##### PARALLEL PROCESSING
    #Leave one or two cores to avoid overload your computer
    cluster <- makeCluster(corenumber)
    registerDoParallel(cluster)

    ##########
    #progress bar
    #########
    registerDoSNOW(cluster)
    iterations <- as.numeric(length(loopseq_z))
    pb <- progress_bar$new(
      format = "letter = :letter [:bar] :elapsed | eta: :eta",
      total = iterations,
      width = 60)

    progress_letter <- rep(LETTERS[1:10], 1000)  # token reported in progress bar

    # allowing progress bar to be used in foreach -----------------------------
    progress <- function(n){
      pb$tick(tokens = list(letter = progress_letter[n]))
    }

    opts <- list(progress = progress)

    gc()

    message("loop starts now")

    system.time(
      a <- foreach (z = 1:length(loopseq_z), .options.snow = opts) %dopar%{
        zarray <- array(0, c(length(loopseq_i), length(loopseq_j)))
        for (l in 1:as.numeric(length(loopseq_i))) {
          for (j in 1:as.numeric(length(loopseq_j))) {
            #get region of interest
            i_begin <- loopseq_i[l] - floor(voisize/2)
            i_end <- loopseq_i[l] + floor(voisize/2)
            j_begin <- loopseq_j[j] - floor(voisize/2)
            j_end <- loopseq_j[j] + floor(voisize/2)
            z_begin <- loopseq_z[z] - floor(voisize/2)
            z_end <- loopseq_z[z] + floor(voisize/2)
            tempvoi <- as.array(big_array[i_begin:i_end, j_begin:j_end,as.numeric(z_begin):as.numeric(z_end)]) #starts at top left corner of matrix

            #check if the pixel at the centre of the box is 0
            voi_mid <- ceiling(0.5*voisize)
            if (tempvoi[voi_mid, voi_mid, voi_mid] == 0) {
              #if it is 0 then the pixel falls outside the bone and we dont want to calculate anything
              # if 0 -> assign 0 and move on to next voi
              zarray[l, j] <- -0.001
            } else {
              #make the voi spherical instead of cubic
              tempvoi[,,] <- tempvoi[,,] * empty_cube[,,]

              marrowcount <- length(tempvoi[tempvoi==1]) #count number of 1's (marrow)
              bonecount <- length(tempvoi[tempvoi==2]) #count number of 2's (bone)
              tempbvtv <- bonecount / (bonecount+marrowcount) #BV/TV in square VOI = bone / (bone+marrow)
              zarray[l, j] <- tempbvtv
            }
            gc()
          }
        }
        bvtv_map[,, z] <- bvtv_map[,, z] + zarray[,]
      }
    )

    gc()

    #mat_array <- do.call(rbind, a)
    ra<-simplify2array(a) #collection of matrices in an array
    dims_ra <- dim(ra)

    #Stop cluster
    stopCluster(cluster) #stop parallel processing


    #create directory to export the bvtv maps
    dir.create(paste(expfolder,"//BVTV_maps", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

    message("exporting output")

    #export views of the BV/TV map in three planes through the center of the image stack
    p <- hcl.colors(50, "Spectral", rev = TRUE) #colours for images

    png(file= paste(expfolder,"//BVTV_maps//" ,folderlist[i],"_BVTV_map_xy_voi_diam", voi_diameter_mm,"mm_interval", voi_interval_mm,"mm.png", sep=""), width=350, height=350)
    image(ra[,,0.5*dims_ra[3]],col=p, main = paste(folderlist[i], " xy BV/TV, marrow", sep=""))
    dev.off()

    png(file= paste(expfolder,"//BVTV_maps//" ,folderlist[i],"_BVTV_map_xz_voi_diam", voi_diameter_mm,"mm_interval", voi_interval_mm,"mm.png", sep=""), width=350, height=650)
    image(ra[,0.5*dims_ra[2],],col=p, main = paste(folderlist[i], " xy BV/TV, marrow", sep=""))
    dev.off()

    png(file= paste(expfolder,"//BVTV_maps//" ,folderlist[i],"_BVTV_map_yz_voi_diam", voi_diameter_mm,"mm_interval", voi_interval_mm,"mm.png", sep=""), width=350, height=650)
    image(ra[0.5*dims_ra[1],,],col=p, main = paste(folderlist[i], " xy BV/TV, marrow", sep=""))
    dev.off()


    #save as nifti
    print(paste("Exporting as nifti - joptool bvtv of : ", folderlist[i]))
    dir.create(paste(expfolder,"//BVTV_maps//" ,folderlist[i],"_BVTV", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    trabnifti<-as.nifti(ra)
    trabnifti@pixdim[2:4]<-c(voxelsize,voxelsize,voxelsize)
    writeNIfTI(nim=trabnifti,filename=paste(expfolder,"//BVTV_maps//",folderlist[i],"_BVTV_voi_diam", voi_diameter_mm,"mm_interval", voi_interval_mm,"mm", sep=""),gzipped=F)
    gc()
    closeAllConnections()


    ########
    ## pointclouds
    ########

    tb_array <- trabnifti
    dims <- dim(tb_array)

    df_out <- data.frame(matrix(nrow=dims[1]*dims[2]*dims[3],ncol=4))
    colnames(df_out) <- c("x","y","z","BVTV")


    output_row_number <- dims[1]*dims[2]*dims[3]


    out_list <- list()

    row_counter <- 1

    for (h in 1:dims[1]) {
      #loop through rows
      for (j in 1:dims[2]) {
        #loop through columns
        for (k in 1:dims[3]) {
          #loop through slices
          df_out <- data.frame(matrix(nrow=1,ncol=4))
          colnames(df_out) <- c("x","y","z","BVTV")
          df_out$x <- h * voi_interval_mm
          df_out$y <- j * voi_interval_mm
          df_out$z <- k * voi_interval_mm
          df_out$BVTV <- tb_array[h,j,k]

          out_list[[row_counter]] <- df_out

          row_counter <- row_counter + 1
        }
      }
    }


    d <- ldply(out_list, data.frame)

    d_cleaned <- d %>% filter(d$BVTV > -0.0001) #remove outside pixels

    d_cleaned$x <- d_cleaned$x - voi_interval_mm
    d_cleaned$y <- d_cleaned$y - voi_interval_mm
    d_cleaned$z <- d_cleaned$z - voi_interval_mm

    d_cleaned$rBVTV <- d_cleaned$BVTV / mean(d_cleaned$BVTV) # calculate rBVTV (Dunmore et al. 2019)

    #range bvtv between minimum (0) and 95% maximum (1) values
    d_cleaned$BVTV_minmax <- d_cleaned$BVTV - min(d_cleaned$BVTV)
    bvtv95 <- quantile(d_cleaned$BVTV_minmax, probs = c(0.95))
    d_cleaned$BVTV_minmax <- d_cleaned$BVTV_minmax / bvtv95



    #export results as CSV
    write.csv(d_cleaned, file = paste(expfolder,"//",folderlist[i],"_BVTV_voi_diam", voi_diameter_mm,"mm_interval", voi_interval_mm,"mm","_pointcloud.csv", sep=""), row.names = FALSE )
  }


}


