#' Get and Process GBIF Data
#' @param data_request GBIF download key
#' @importFrom rgbif occ_download_meta
#' @importFrom bigreadr nlines big_fread1
#' @importFrom dplyr select distinct
GetGBIFData <- function(data_request) {
  # GBIF meta returns a key; ensure we use the character string
  key <- as.character(data_request)
  zip_file <- paste0(key, ".zip")
  
  if (file.exists(zip_file)) {
    # Extract files and identify the data file (usually occurrence.txt)
    extracted_files <- unzip(zip_file, list = TRUE)$Name
    data_file <- extracted_files[grep("\\.(csv|txt)$", extracted_files)][1]
    
    unzip(zip_file, files = data_file)
    
    message("Processing ", data_file, "...")
    
    occs <- bigreadr::big_fread1(
      data_file,
      every_nlines = 1000000,
      integer64 = "double",
      .transform = function(x) {
        x %>%
          dplyr::select(
            gbifID, family, species, decimalLongitude, decimalLatitude,
            elevation, eventDate, countryCode, coordinateUncertaintyInMeters
          ) %>%
          dplyr::distinct(species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)
      }
    )
    
    # Clean up extracted file and zip
    file.remove(data_file)
    file.remove(zip_file) 
    
    # Final global
    occs <- dplyr::distinct(occs, species, decimalLongitude, decimalLatitude, eventDate, .keep_all = TRUE)
    return(occs)
  } else {
    stop("Zip file not found. Ensure occ_download_get finished successfully.")
  }
}

DownloadGBIF <- function(key, user, user.email, pwd, custom.shp) {
  if (is.null(custom.shp)) stop("A valid spatial object (sf/Spatial) is required.")
  
  # 1. Convert shape to WKT
  aoi_wkt <- rgbif::gbif_bbox2wkt(bbox = sf::st_bbox(custom.shp))
  
  # 2. Trigger download - CHANGE FORMAT TO "DWCA"
  req <- rgbif::occ_download(
    rgbif::pred_in("taxonKey", key),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred_within(aoi_wkt),
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
  if(!dir.exists("data/dwc_raw")) dir.create("data/dwc_raw", recursive = TRUE)
  
  download_info <- rgbif::occ_download_get(req, path = "data/dwc_raw", overwrite = TRUE)
  
  # DwC-A contains multiple files, so we unzip into a folder named by the Key
  unzip_dir <- file.path("data/dwc_raw", as.character(req))
  unzip(download_info, exdir = unzip_dir)
  
  # In DwC-A, the main data is ALWAYS 'occurrence.txt'
  occurrence_file <- file.path(unzip_dir, "occurrence.txt")
  
  return(list(
    key = as.character(req),
    file_path = occurrence_file,
    all_files = list.files(unzip_dir, full.names = TRUE),
    citation = rgbif::gbif_citation(as.character(req))
  ))
}

