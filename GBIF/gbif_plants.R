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
sf::sf_use_s2(FALSE)
# Define hemispheres
west <- st_as_sfc(st_bbox(c(xmin = -180, ymin = -90, xmax = 0, ymax = 90), crs = 4326))
east <- st_as_sfc(st_bbox(c(xmin = 0, ymin = -90, xmax = 180, ymax = 90), crs = 4326))
# download mountains areas
mountains <- st_read("GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp") %>%
  st_make_valid()
# Split problematic polygons by hemisphere
prob <- mountains[c(1, 12, 100), ]
prob_west <- st_intersection(prob, west)
prob_east <- st_intersection(prob, east)
# Recombine with clean polygons
mountains_clean <- mountains[-c(1, 12, 100), ]
mountains_buf <- rbind(mountains_clean, prob_west, prob_east) %>%
  st_convex_hull() %>%
  st_union() %>%
  st_make_valid() %>%
  wk::wk_orient()

# Get Taxon Key
taxon_key <- name_backbone(name = "Aves")$classKey

# Note: Ensure it returns a dataframe or a path to the CSV
gbif_res <- DownloadGBIF(taxon_key, gbif_user, gbif_email, gbif_pwd, mountains_buf)
file_path <- gbif_res$file_path # extract csv path

# 2. Convert GBIF CSV to Parquet (The "Arrow" Way) ----
# Do this once! It will shrink a 10GB CSV into a ~1GB Parquet file.
parquet_dir <- "GBIF/data/Aves_parquet"
GBIF_to_Parquet(file_path = file_path, parquet_dir = parquet_dir)

# 3. Data cleaning
ds <- arrow::open_dataset("GBIF/data/Aves_parquet")
clean_path <- "GBIF/data/Aves_parquetclean"
years <- 2009:2026

if (!dir.exists(clean_path)) dir.create(clean_path, recursive = TRUE)

for (y in years) {
  
  cat("Processing year:", y, "\n")
  
  # 3.1. Download year per year
  occ <- ds %>%
    filter(year == y) %>%
    
    # Filtering on basisOfRecord isn't necessary because it's done earlier (DownloadGBIF)
    filter(taxonRank %in% c("SPECIES", "SUBSPECIES")) %>%
    filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 250) %>%
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
  message("Processing to spatial filtering.")
  # convert data frame to sf object
  occ_sf <- st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
  # join our data frame with mountains polygons (keep occurrences only within mountains polygons)
  occ_sf <- st_join(occ_sf, mountains, join = st_intersects, left = FALSE)
  
  if (nrow(occ_sf) == 0) next
  
  # back to dataframe
  occ <- st_drop_geometry(occ_sf)
  
  message("Writing to Parquet.")
  # 3.4. re-save in parquet again
  arrow::write_parquet(occ, file.path(clean_path, paste0("occ_", y, ".parquet")))
  
  rm(occ)
  rm(occ_sf)
  rm(occ_cen)
  gc()
  
}

