#reads a 3D heatmap as a nifti file and exports a pointcloud as csv file with and without zeros
trabmap_to_pointcloud <- function(folder = "", writefolder = "", spacing_mm = 2.5 ) {

  folderlist <- trabmap_folder_list(folder)
  for (i in 1:length(folderlist)) {
    tb_array <- readNIfTI(paste(folder, folderlist[i], sep=""))


    dims <- dim(tb_array)

    df_out <- data.frame(matrix(nrow=dims[1]*dims[2]*dims[3],ncol=4))
    colnames(df_out) <- c("x","y","z","BVTV")


    output_row_number <- dims[1]*dims[2]*dims[3]


    out_list <- list()

    row_counter <- 1

    for (i in 1:dims[1]) {
      #loop through rows
      for (j in 1:dims[2]) {
        #loop through columns
        for (k in 1:dims[3]) {
          #loop through slices
          df_out <- data.frame(matrix(nrow=1,ncol=4))
          colnames(df_out) <- c("x","y","z","BVTV")
          df_out$x <- i * spacing_mm
          df_out$y <- j * spacing_mm
          df_out$z <- k * spacing_mm
          df_out$BVTV <- tb_array[i,j,k]

          out_list[[row_counter]] <- df_out

          row_counter <- row_counter + 1
        }
      }
    }


    d <- ldply(out_list, data.frame)

    expfolder <- c(paste(writefolder,"//", folderlist[i],  sep = ""))
    dir.create(expfolder, showWarnings = TRUE, recursive = FALSE, mode = "0777")


    write.csv(d, file = paste(expfolder,"//",folderlist[i],"_pointcloud.csv", sep=""), row.names = FALSE )

    d_cleaned <- d %>% filter(d$BVTV > 0)


    write.csv(d, file = paste(expfolder,"//",folderlist[i],"_pointcloud_no_zeros.csv", sep=""), row.names = FALSE )

  }

}



