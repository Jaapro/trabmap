#library(remotes)
#remotes::install_github("Jaapro/trabmap")
#library(trabmap)
#trabmap_load()

location <- "C://CT_scans//trabmap_input//"
out <- "C://CT_scans//trabmap_output"

vol <- trabmap_readtiffstack(location, cores=1)

trabmap_fill_trinary(folder = location, writefolder = "C://CT_scans//trabmap_output", voxelsizelist = c(0.03))


trabmap_fill_trinary_chopped(folder = location, writefolder = "C://CT_scans//trabmap_output", voxelsizelist = c(0.03))

trabmap_map(folder = location, writefolder = out, voxelsizelist = c(0.03), voi_diameter_mm = 2, voi_interval_mm = 0.5 )



trabmap_combine_nifti(folder = location, writefolder = out, voxelsize = c(0.3))


trabmap_combine_nifti_b(folder = location, writefolder = out, voxelsize = c(0.3))



##debugging
folderlist <- trabmap_folder_list(location)
i<- 1
tb_array <- readNIfTI(paste(location, folderlist[i], sep=""))
dims <- dim(tb_array)


trabnifti <- as.nifti(tb_array)
trabnifti@data_type<- "int8"
trabnifti@datatype <- "int16"
convert.datatype()

writeNIfTI(nim=trabnifti,filename="C://CT_scans//trabmap_output//test")

tn <- datatyper(img = tb_array, type_string = "INT8", datatype = "INT8")
writeNIfTI(nim=tn,filename="C://CT_scans//trabmap_output//test_b")

writeNIfTI2()





m <- array(0, dim = c(1000, 70, 1, 1000))
format(object.size(m), units = "auto")

m3 <- array(raw(0), dim = c(1000, 70, 1, 1000))
format(object.size(m3), units = "auto")

format(object.size(tb_array), units = "auto")
tb_int <- as.integer(tb_array)
format(object.size(tb_int), units = "auto")
### transforming to raw and then to int halves filesize!

df <- as.array(tb_array)
format(object.size(df), units = "auto")
image(df[,,8])

rawdf <- as.raw(df)
rawarray <- array(rawdf, c(722, 528, 36))
format(object.size(rawarray), units = "auto")
image(rawarray[,,5])

rawint <- as.integer(rawarray)
rawint <- array(rawint, c(722,528,36))
format(object.size(rawint), units = "auto")
image(rawint[,,5])

