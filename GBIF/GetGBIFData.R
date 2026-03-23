# ------
# Function to download GBIF data from gbif.com into an unzipped file
DownloadGBIF <- function(key, user, user.email, pwd, custom.shp) {
  if (is.null(custom.shp)) stop("A valid spatial object (sf/Spatial) is required.")
  
  # 1. Convert shape to WKT
  aoi_wkt <- sf::st_as_text(custom.shp)
  
  # 2. Trigger download - CHANGE FORMAT TO "DWCA"
  req <- rgbif::occ_download(
    # rgbif::pred_in("year", 1980), # Add only if the file is too big
    rgbif::pred_in("taxonKey", key),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    rgbif::pred("occurrenceStatus", "PRESENT"),
    rgbif::pred_within(aoi_wkt),
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
  unzip(download_info, files = "occurrence.txt", exdir = unzip_dir)
  
  # Remove the original zip file
  file.remove(download_info)
  
  # In DwC-A, the main data is ALWAYS 'occurrence.txt'
  occurrence_file <- file.path(unzip_dir, "occurrence.txt")
  
  return(list(
    key = as.character(req),
    file_path = occurrence_file,
    citation = rgbif::gbif_citation(as.character(req))
  ))
}

# ------
# Function to convert GBIF CSV to Parquet (The "Arrow" Way) 
GBIF_to_Parquet <- function(file_path, parquet_dir) {
  
  # 1. Create output directory if needed
  if (!dir.exists(parquet_dir)) {
    dir.create(parquet_dir, recursive = TRUE)
  }
  
  message("Starting chunked read of GBIF CSV...")
  
  chunk_id <- 0
  
  # 2. Read + process in chunks + write to Parquet
  bigreadr::big_fread1(file_path, every_nlines = 1e6, integer64 = "double", 
                               nThread = parallel::detectCores()-1, verbose = TRUE, 
                               .transform = function(x) {
                                 chunk_id <- chunk_id + 1
      message("Processing chunk of ", chunk_id, " of ", nrow(x), " rows...")
      processed <- x %>%
        dplyr::select(
          species, scientificName, kingdom, phylum, class, order, family, genus, 
          taxonRank, occurrenceID, basisOfRecord, eventDate, year, countryCode,
          decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters
        ) %>%
        dplyr::distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)
      
      arrow::write_dataset(processed, path = parquet_dir, format = "parquet", partitioning = "year",
                           basename_template = paste0("chunk_", chunk_id, "_part{i}.parquet"))
      
      return(NULL)
    })
  
  # 4. Cleanup 
  unlink(dirname(file_path), recursive = TRUE)
  gc()
  
  message("Conversion to Parquet complete!")
}
