pad_image_stack <- function(img_stack, pad_thickness = 1) {
  dims <- dim(img_stack)
  # Create a new empty array with extra padding layers
  padded_stack <- array(0, dim = dims + 2 * pad_thickness)
  # Place the original image stack in the center
  padded_stack[(pad_thickness + 1):(dims[1] + pad_thickness),
               (pad_thickness + 1):(dims[2] + pad_thickness),
               (pad_thickness + 1):(dims[3] + pad_thickness)] <- img_stack
  return(padded_stack)
}

# Function to convert binary image stack to a 3D surface mesh and save it
trabmap_stack_to_ply <- function(image_stack, voxelsize, output_ply_path, reduction_factor = 0.25, pad_thickness = 10, resize_factor = 0.5) {
  # Load binary image stack
  # Assumes the image stack is in .tiff format and the images are binary (0 and 1)
  voxel_dims <- c(voxelsize, voxelsize, voxelsize)

  img_stack <- image_stack@.Data[seq(1, dim(image_stack)[1], by = 1/resize_factor),
                                 seq(1, dim(image_stack)[2], by = 1/resize_factor),
                                 seq(1, dim(image_stack)[3], by = 1/resize_factor)]

  # Replace with your original voxel size
  new_voxel_size <- voxel_dims / resize_factor

  img_stack <- nifti(array(img_stack, dim = dim(img_stack)),
                     pixdim = c(1, new_voxel_size))

  # Ensure the image stack is binary (0 and 1)
  img_stack[img_stack > 0] <- 1

  img_stack <- pad_image_stack(img_stack, pad_thickness)

  # Create 3D surface mesh
  mesh <- vcgIsosurface(img_stack, threshold = 0.5)
  rm(img_stack)

  #set voxelsize
  mesh$vb[1:3, ] <- mesh$vb[1:3, ] * new_voxel_size

  # Simplify the mesh
  simplified_mesh <- vcgQEdecim(mesh, percent = reduction_factor)
  rm(mesh)

  # Create a flip matrix for X and Y axes
  simplified_mesh$vb[1, ] <- -simplified_mesh$vb[1, ] #flip x axis
  simplified_mesh$vb[2, ] <- -simplified_mesh$vb[2, ] #flip x axis

  # Attempt to close holes in the mesh
  closed_mesh <- vcgClean(simplified_mesh, sel=c(2,3,4), iterate = TRUE)
  rm(simplified_mesh)

  # Export the mesh to ASCII .ply format
  vcgPlyWrite(closed_mesh, output_ply_path, ascii = TRUE)

  rm(closed_mesh)

  cat("Mesh has been successfully saved to", output_ply_path, "\n")
}

