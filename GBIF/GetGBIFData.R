# ------
# Function to download GBIF data from gbif.com into an unzipped file
DownloadGBIF <- function(key, user, user.email, pwd, custom.shp) {
  if (is.null(custom.shp)) stop("A valid spatial object (sf/Spatial) is required.")
  
  # 1. Convert shape to WKT
  aoi_wkt <- rgbif::gbif_bbox2wkt(bbox = sf::st_bbox(custom.shp))
  
  # 2. Trigger download - CHANGE FORMAT TO "DWCA"
  req <- rgbif::occ_download(
    rgbif::pred_in("taxonKey", key),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    rgbif::pred_within(aoi_wkt),
    rgbif::pred_in("taxonRank", c("SPECIES", "SUBSPECIES")),
    rgbif::pred_in("basisOfRecord",
                   c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN")),
    format = "DWCA", # <--- This is the Darwin Core Archive format
    user = user, pwd = pwd, email = user.email
  )
  
  # 3. Polling for completion
  repeat {
    status <- rgbif::occ_download_meta(req)$status
    message("Status: ", status, " (Key: ", as.character(req), ")")
    if (status %in% c("SUCCEEDED", "KILLED", "CANCELLED")) break
    Sys.sleep(45)
  }
  
  if (status != "SUCCEEDED") stop("GBIF Download failed.")
  
  # 4. Retrieval and Unzipping
  if(!dir.exists("GBIF/data/dwc_raw")) dir.create("GBIF/data/dwc_raw", recursive = TRUE)
  
  download_info <- rgbif::occ_download_get(req, path = "GBIF/data/dwc_raw", overwrite = TRUE)
  
  # DwC-A contains multiple files, so we unzip into a folder named by the Key
  unzip_dir <- file.path("GBIF/data/dwc_raw", as.character(req))
  unzip(download_info, exdir = unzip_dir)
  
  # Remove the original zip file
  file.remove(download_info)
  
  # In DwC-A, the main data is ALWAYS 'occurrence.txt'
  occurrence_file <- file.path(unzip_dir, "occurrence.txt")
  
  return(list(
    key = as.character(req),
    file_path = occurrence_file,
    all_files = list.files(unzip_dir, full.names = TRUE),
    citation = rgbif::gbif_citation(as.character(req))
  ))
}

# ------
# Function to convert GBIF CSV to Parquet (The "Arrow" Way) 
GBIF_to_Parquet <- function(file_path, parquet_dir) {
  
  # 1. Create output directory if needed
  if (!dir.exists(parquet_dir)) {
    
    dir.create(parquet_dir, recursive = TRUE)
  
  message("Starting chunked read of GBIF CSV...")
  
  # 2. Read + process in chunks
  occs <- bigreadr::big_fread1(file_path, every_nlines = 1e6, integer64 = "double",
    verbose = TRUE, .transform = function(x) {
      message("Processing chunk of ", nrow(x), " rows...")
      x %>%
        dplyr::select(
          species, scientificName, kingdom, phylum, class, order, family, genus, 
          taxonRank, occurrenceID, basisOfRecord, eventDate, year, countryCode,
          decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters
        ) %>%
        dplyr::distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)
    }
  )
  
  # 3. Write to Parquet
  message("Writing to Parquet...")
  arrow::write_dataset(occs, path = parquet_dir, format = "parquet", partitioning = "year")
  
  # 4. Cleanup 
  unlink(dirname(file_path), recursive = TRUE)
  rm(occs)
  gc()
  
  message("Conversion to Parquet complete!")
  
  }
}
