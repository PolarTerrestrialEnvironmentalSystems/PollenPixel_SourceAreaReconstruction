#script to download input data
rm(list = ls())

# unzip the input dir from zenodo
data_exists <- dir.exists("input")
if(!data_exists){
  options(timeout = 1e6)
  
  #download from Zenodo
  download.file("https://zenodo.org/records/21129432/files/input.tar.gz?download=1", 
                destfile = "input.tar.gz",
                mode = "wb",
                method = "libcurl")
  
  #unpack
  untar("input.tar.gz")
  
  #remove tar.gz file
  unlink("input.tar.gz",
         recursive = TRUE)
}

# download CHELSA bioclim data (wget)
# check whether CHELSA data has already been downloaded
clim_files <- list.files(path = "input/data/global_maps",
                         pattern = "CHELSA")

if(length(clim_files) != 22){
  file_links <- readLines("input/data/global_maps/clim_data.txt")
  for(file in file_links){
    file_name <- unlist(strsplit(file, "/"))[10]
    if(file.exists(file.path("input/data/global_maps",file_name))) next
    options(timeout = 1e6)
    download.file(file,
                  destfile = file.path("input/data/global_maps",
                                       file_name),
                  mode = "wb",
                  method = "libcurl")
  }
}

# download ETOPO topography (wget)
topo_there <- file.exists("input/data/global_maps/ETOPO1_Ice_g_geotiff.tif")
if(!topo_there){
  topo_link <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/georeferenced_tiff/ETOPO1_Ice_g_geotiff.zip"
  #download
  download.file(topo_link,
                destfile = "input/data/global_maps/ETOPO1_Ice_g_geotiff.zip",
                method = "libcurl",
                mode = "wb")
  #unzip
  unzip("input/data/global_maps/ETOPO1_Ice_g_geotiff.zip",
        "ETOPO1_Ice_g_geotiff.tif",
        exdir = "input/data/global_maps")
  
  #delete zipped version
  unlink("input/data/global_maps/ETOPO1_Ice_g_geotiff.zip",
         recursive = TRUE)
}

# check and/or install necessary packages
packages <- readLines("scripts/libraries.txt")
check_and_load <- function(package_name){
  is_there <- require(package_name,
                      character.only = TRUE)
  if(!is_there){
    print("installing ",package_name)
    install.packages(package_name)
    require(package_name,
            character.only = TRUE)
  }
}

sapply(packages,
       check_and_load,
       USE.NAMES = FALSE)
sf_use_s2(FALSE)
