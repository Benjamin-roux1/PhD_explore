# GBIF DATA
# 17-03-2026
# Benjamin Roux

# 1. Libraries ----
library(tidyverse)
library(arrow)      
library(rgbif)
library(sf)
library(CoordinateCleaner)
library(data.table)
library(dplyr)
library(bigreadr)
source("GBIF/GetGBIFData.R")

# Increase memory overhead for large spatial joins if needed
options(future.globals.maxSize = 8000 * 1024^2) 

# 2. GBIF Authentication ----
gbif_user  <- Sys.getenv("GBIF_USER") %||% rstudioapi::askForPassword("GBIF Username")
gbif_email <- Sys.getenv("GBIF_EMAIL") %||% rstudioapi::askForPassword("GBIF Email")
gbif_pwd   <- Sys.getenv("GBIF_PWD") %||% rstudioapi::askForPassword("GBIF Password")

# Load Mountain Shapefile
mountains <- st_read("GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp") %>%
  st_make_valid() #

# Get Taxon Key
taxon_key <- name_backbone(name = "Aves", rank="class")$classKey

# Note: Ensure it returns a dataframe or a path to the CSV
gbif_res <- DownloadGBIF(taxon_key, gbif_user, gbif_email, gbif_pwd, mountains)
file_path <- gbif_res$file_path # extract csv path

# 2. Convert GBIF CSV to Parquet (The "Arrow" Way) ----
# Do this once! It will shrink a 10GB CSV into a ~1GB Parquet file.
parquet_dir <- "GBIF/data/gbif_parquet"
GBIF_to_Parquet(file_path = file_path, parquet_dir = parquet_dir)

# 3. Data cleaning
ds <- arrow::open_dataset("GBIF/data/gbif_parquet")

clean_path <- "GBIF/data/gbif_parquetclean"
years <- 1980:2026

if (!dir.exists(clean_path)) dir.create(clean_path, recursive = TRUE)

for (y in years) {
  
  cat("Processing year:", y, "\n")
  
  # 3.1. Download year per year
  occ <- ds %>%
    filter(year == y) %>%

    filter(taxonRank %in% c("SPECIES", "SUBSPECIES"), basisOfRecord %in% c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN")) %>%
    collect()

  if (nrow(occ) == 0) next
  
  # 3.2. CoordinateCleaner
  occ$val <- cc_val(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$zero <- cc_zero(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$equ <- cc_equ(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$cap <- cc_cap(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$inst <- cc_inst(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ_cen <- cc_cen(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ <- occ[occ$val & occ$zero & occ$equ & occ$cap & occ$inst & occ_cen, ]
  occ <- cc_dupl(occ, lon = "decimalLongitude", lat = "decimalLatitude")
  
  # 3.3. Spatial filtering: keep occurrences only inside the mountains ranges
  # convert data frame to sf object
  occ_sf <- st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
  # join our data frame with mountains polygons (keep occurrences only within mountains polygons)
  occ_sf <- st_join(occ_sf, mountains, join = st_intersects, left = FALSE)
  
  if (nrow(occ_sf) == 0) next
  
  # back to dataframe
  occ <- st_drop_geometry(occ_sf)
  
  # 3.4. re-save in parquet again
  arrow::write_parquet(occ, file.path(clean_path, paste0("occ_", y, ".parquet")))
  
  rm(occ)
  rm(occ_sf)
  gc()
}

