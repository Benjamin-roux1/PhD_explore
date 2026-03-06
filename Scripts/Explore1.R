install.packages("terra")
install.packages("sf")
install.packages("leaflet")
install.packages("mapview")

library(terra); library(sf); library(leaflet); library(mapview)

GMBA <- st_read("GMBA_inventory_v2.0_standard_300/GMBA_inventory_v2.0_standard_300.shp")

mapview(GMBA, zcol = "MapName", legend = FALSE)

GMBA_vect <- terra::vect(GMBA)
plot(GMBA_vect) 

CHELSA_clim <- read.table("CHELSA/CHELSA_clim.txt")[,1]
NAME <- unlist(lapply(strsplit(CHELSA_clim,"/"),function(x) x[10]))

for (i in 1:length(CHELSA_clim)){
  print(paste(i,length(CHELSA_clim),sep=" out of "))
  shell(cmd=paste('wget -nc -c -nd ',CHELSA_clim[i],sep=""))
  tmp <- terra::mask(terra::rast(NAME[i]),GMBA_vect)
  terra::writeRaster(tmp,file.path("CHELSA", paste0("GMBA_", NAME[i])),overwrite=T)
  file.remove(NAME[i])
}

tif_files <- list.files("CHELSA/", pattern = "\\.tif$", full.names = TRUE)
rasters_stack <- terra::rast(tif_files)
plot(rasters_stack)
