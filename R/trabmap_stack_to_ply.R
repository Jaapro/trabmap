install.packages("EBImage")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EBImage")

install.packages("Rvcg")
install.packages("rgl")


library(EBImage)
library(Rvcg)
library(rgl)
library(oro.nifti)

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
convert_and_simplify_mesh <- function(image_stack_path, voxel_dims, output_ply_path, reduction_factor = 0.05) {
  # Load binary image stack
  # Assumes the image stack is in .tiff format and the images are binary (0 and 1)
  img_stack <- readNIfTI(image_stack_path)

  # Ensure the image stack is binary (0 and 1)
  img_stack[img_stack > 0] <- 1

  img_stack <- pad_image_stack(img_stack, pad_thickness)

  # Create 3D surface mesh
  mesh <- vcgIsosurface(img_stack, threshold = 0.5)

  #set voxel size
  mesh$vb[1:3, ] <- mesh$vb[1:3, ] * voxel_dims


  # Simplify the mesh
  simplified_mesh <- vcgQEdecim(mesh, percent = reduction_factor)

  # Attempt to close holes in the mesh
  closed_mesh <- vcgClean(simplified_mesh, sel=c(2,3,4), iterate = TRUE)

  # Export the mesh to ASCII .ply format
  vcgPlyWrite(simplified_mesh, output_ply_path, ascii = TRUE)

  cat("Mesh has been successfully saved to", output_ply_path, "\n")
}

# Parameters
image_stack_path <- "D://CT_scans//trabmap//input//AIM11444_pongo_pygmaeus_HUM_dist_FILL_trinary.nii"  # Path to binary image stack
voxel_dims <- c(0.5, 0.5, 0.5)  # Set your voxel dimensions (x, y, z) in mm or other units
output_ply_path <- "D://CT_scans//trabmap//input//AIM11444_pongo_pygma_output_mesh.ply"  # Output .ply file path
reduction_factor <- 0.05  # Adjust reduction factor (0.5 reduces the mesh to 50% of original)
pad_thickness <- 10  # Number of voxel layers to pad the image stack (increase if needed)

# Run the function
# Run the function
convert_and_simplify_mesh(image_stack_path, voxel_dims, output_ply_path, reduction_factor)
