CHELSA_clim <- read.table("CHELSA/CHELSA_clim.txt")[,1]
NAME <- unlist(lapply(strsplit(CHELSA_clim,"/"),function(x) x[10]))

for (i in 1:length(CHELSA_clim)){
  print(paste(i,length(CHELSA_clim),sep=" out of "))
  shell(cmd=paste('wget -nc -c -nd ',CHELSA_clim[i],sep=""))
  tmp <- terra::mask(terra::rast(NAME[i]),GMBA_vect)
  terra::writeRaster(tmp,file.path("CHELSA", paste0("GMBA_", NAME[i])),overwrite=T)
  file.remove(NAME[i])
}
