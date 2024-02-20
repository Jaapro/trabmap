trabmap_fill <- function(mask, steps = 4, strel = 13){
  require(EBImage)

  d_e_strel <-  makeBrush(size=strel,shape="disc")
  d_e_steps <- steps

  if(length(d_e_steps)==0){
    mask<-EBImage::fillHull(mask)
  } else {
    for(j in 1:d_e_steps){mask<-EBImage::dilate(mask,d_e_strel)}
    #mask<-EBImage::closing(mask,d_e_strel)
    mask<-EBImage::fillHull(mask)
    for(j in 1:d_e_steps){mask<-EBImage::erode(mask,d_e_strel)}
  }
  return(mask)
}
