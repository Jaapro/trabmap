# make a 3d sphere to use as a structuring element for erosion etc

trabmap_strel3d <- function(diameter) {
  empty_cube <- array(0,c(diameter,diameter,diameter))
  loop_length <- ceiling(diameter/2)
  expand_factor <- 0
  for (l in 1:loop_length){
    if (l != loop_length){
      empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor) , l] <- empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor) , l] + makeBrush(size=l+expand_factor, shape="disc") #start of stack
      empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor)  , diameter - expand_factor] <- empty_cube[(loop_length-expand_factor):(loop_length+expand_factor) ,(loop_length-expand_factor):(loop_length+expand_factor) , diameter - expand_factor] + makeBrush(size=l+expand_factor, shape="disc") #end of stack
      expand_factor <- expand_factor + 1
    }
    else {
      empty_cube[ , , loop_length] <- empty_cube[ , , loop_length] + makeBrush(size=diameter, shape="disc") #start of stack
      }
  }
  return(empty_cube)
}
