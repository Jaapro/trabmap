#' Read TIFF image stack
#'
#' FUNCTION COPIED FROM ALESSIO VENEZIANO'S INDIANABONES PACKAGE
#'
#' Reads TIFF image stacks into an array. It uses the function \code{readTIFF} in the package 'tiff' and runs in parallel using the packages 'foreach' and 'doParallel'.
#'
#' @param path path of the folder where the TIFF images are stored
#' @param cores number of cores for parallel computation (default=2)
#'
#' @return An array with dimensions depending on the number of channels in the images. The first two dimensions indicate the image width and height in pixels, and the last dimension is the number of images in the stack.
#'
#' @author Alessio Veneziano
#'
#' @export

trabmap_readtiffstack<-function(path,cores = 1){
  require(tiff)
  require(doParallel)
  require(foreach)

  files<-list.files(path)
  files<-paste(path,files,sep="/")

  registerDoParallel(cores=cores)
  suppressWarnings(
    im<-foreach(i=files) %dopar% tiff::readTIFF(i)
  )
  dims<-dim(im[[1]])
  im<-array(unlist(im),c(dims,length(files)))

  if(length(dims)>2){
    perm <- 1:length(dims)
    perm[1:2]<-perm[2:1]
    perm<-c(perm,4)
    im<-aperm(im,perm)
  } else if( length(dims)==2){
    perm<-c(2,1,3)
    im<-aperm(im,perm)
  }

  return(im)
}
