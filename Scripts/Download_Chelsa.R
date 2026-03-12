# CHELSA download
# Benjamin Roux
# 12-03-2026

# Download CHELSA text file on CHELSA database
# TF <- read.table("CHELSA/CHELSA_test.txt")[,1] #text file from CHELSA database with each files path
# method: must be one of "mask", "crop", "raw".
# vect: SpatVector to mask the CHELSA file (e.g. GMBA mountains ranges)
# ext: manually defined extent --> MYEXT <- terra::ext(-12, 50, 28, 55) # east and north extension
# outdir: output directory

# Function to download CHELSA files according to 3 scenarios: (i) raw file, (ii) cropped or (iii) masked
download_chelsa <- function(TF, method, vect, ext, outdir) {
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # For each file in TF
  for (i in seq_along(TF)) {
    
    NAME <- basename(TF[i])
    print(paste(i, "out of", length(TF)))
    
    # Download with wget
    system(paste("wget -nc -c -nd", TF[i]))
    
    r <- terra::rast(NAME)
    
    if (method == "mask") {
      if (is.null(vect)) stop("vect must be provided for mask method")
      r <- terra::mask(r, vect)
    } else if (method == "crop") {
      if (is.null(ext)) stop("ext must be provided for crop method")
      r <- terra::crop(r, ext)
    } else if (method == "raw") {
      # Do nothing
    } else {
      stop("method must be one of 'mask', 'crop', 'raw'")
    }
    
    # Name the file adding GMBA, should be personalized
    outname <- file.path(outdir, paste0("GMBA_", NAME))
    terra::writeRaster(r, outname, overwrite = TRUE)
    
    file.remove(NAME)
  }
  
  print("All files processed!")
}
