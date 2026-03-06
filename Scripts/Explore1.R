setwd("C:/Users/berou1714/OneDrive - Norwegian University of Life Sciences/Desktop/PhD-explore")

install.packages("terra")
install.packages("sf")
install.packages("leaflet")
install.packages("mapview")

library(terra); library(sf); library(leaflet); library(mapview)

GMBA <- st_read("GMBA_inventory_v2.0_standard_300/GMBA_inventory_v2.0_standard_300.shp")

mapview(GMBA, zcol = "MapName", legend = FALSE)

GMBA_vect <- terra::vect(GMBA)
plot(GMBA_vect) #ej fais une modification ici 

CHELSA_clim <- read.table("CHELSA/CHELSA_clim.txt")[,1]
NAME <- unlist(lapply(strsplit(CHELSA_clim,"/"),function(x) x[10]))

for (i in 1:length(CHELSA_clim)){
  print(paste(i,length(CHELSA_clim),sep=" out of "))
  shell(cmd=paste('wget -nc -c -nd ',CHELSA_clim[i],sep=""))
  tmp <- terra::mask(terra::rast(NAME[i]),GMBA_vect)
  terra::writeRaster(tmp_mask,file.path("CHELSA", paste0("GMBA_", NAME[i])),overwrite=T)
  file.remove(NAME[i])
}

#je m'apelle
#oui oui
#coebnjfdcneb