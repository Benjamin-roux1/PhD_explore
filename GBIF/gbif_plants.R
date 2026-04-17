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
gbif_res <- Download.GBIF(taxon_key, gbif_user, gbif_email, gbif_pwd, mountains_buf)
file_path <- gbif_res$file_path # extract csv path

# 2. Convert GBIF CSV to Parquet (The "Arrow" Way) ----
# Do this once! It will shrink a 10GB CSV into a ~1GB Parquet file.
parquet_dir <- "GBIF/data/Aves_parquet"
GBIF.to.Parquet(file_path = file_path, parquet_dir = parquet_dir)

# 3. Data cleaning
ds <- arrow::open_dataset("GBIF/data/Aves_parquet")
clean_path <- "GBIF/data/Aves_parquetclean"

if (!dir.exists(clean_path)) dir.create(clean_path, recursive = TRUE)
years <- 1980:2026

for (y in years) {
  
  cat("Processing year:", y, "\n")
  
  # 3.1. Download year per year
  occ <- ds %>%
    filter(year == y,  # Filtering on basisOfRecord isn't necessary because it's done earlier (DownloadGBIF)
           taxonRank %in% c("SPECIES", "SUBSPECIES"),
           !is.na(decimalLatitude),
           !is.na(decimalLongitude),
           !is.na(species)) %>%  
    collect()

  if (nrow(occ) == 0) { message("No records for ", y, ", skipping."); next }
  
  # 3.2. CoordinateCleaner
  occ$val <- cc_val(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$zero <- cc_zero(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$equ <- cc_equ(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$cap <- cc_cap(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$inst <- cc_inst(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
  occ$cen <- cc_cen(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
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
  gc()
  
}


# --------
# -------- Same but for very big datasets (process by chunks)
if (!dir.exists(clean_path)) dir.create(clean_path, recursive = TRUE)
years <- 1980:2026

for (y in years) {
  cat("\n── Processing year:", y, "──\n")
  
  # Get all parquet chunks for this year
  year_dir <- file.path("GBIF/data/Aves_parquet", paste0("year=", y))
  if (!dir.exists(year_dir)) { message("No directory for ", y, ", skipping."); next }
  
  chunks <- list.files(year_dir, pattern = "\\.parquet$", full.names = TRUE)
  if (length(chunks) == 0) { message("No chunks for ", y, ", skipping."); next }
  cat("  Found", length(chunks), "chunks\n")
  
  year_results <- vector("list", length(chunks))  # pre-allocate list
  
  for (j in seq_along(chunks)) {
    cat("  Chunk", j, "/", length(chunks), "\n")
    
    # ── 3.1. Read single chunk ──────────────────────────────────────────────
    occ <- arrow::read_parquet(chunks[j]) %>%
      filter(taxonRank %in% c("SPECIES", "SUBSPECIES"),
             !is.na(decimalLatitude),
             !is.na(decimalLongitude),
             !is.na(species))
    
    if (nrow(occ) == 0) next
    
    # ── 3.2. CoordinateCleaner ──────────────────────────────────────────────
    flags <- data.frame(
      val  = cc_val(occ,  lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged"),
      zero = cc_zero(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged"),
      equ  = cc_equ(occ,  lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged"),
      cap  = cc_cap(occ,  lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged"),
      inst = cc_inst(occ, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged"),
      cen  = cc_cen(occ,  lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")
    )
    occ <- occ[rowSums(!flags) == 0, ]
    occ <- cc_dupl(occ, lon = "decimalLongitude", lat = "decimalLatitude")
    
    if (nrow(occ) == 0) next
    
    # ── 3.3. Spatial filtering ──────────────────────────────────────────────
    occ_sf <- st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"),
                       crs = 4326, remove = FALSE)
    occ_sf <- st_join(occ_sf, mountains, join = st_intersects, left = FALSE)
    
    if (nrow(occ_sf) == 0) next
    
    year_results[[j]] <- st_drop_geometry(occ_sf)
    rm(occ, occ_sf, flags); gc()
  }
  
  # ── 3.4. Combine all chunks for this year and write once ─────────────────
  year_results <- Filter(Negate(is.null), year_results)  # drop empty chunks
  if (length(year_results) == 0) { message("No clean records for ", y, ", skipping."); next }
  
  occ_year <- data.table::rbindlist(year_results, fill = TRUE)
  
  # Re-run cc_dupl across the full year to catch cross-chunk duplicates
  occ_year <- cc_dupl(occ_year, lon = "decimalLongitude", lat = "decimalLatitude")
  
  cat("  Final records for", y, ":", nrow(occ_year), "\n")
  arrow::write_parquet(occ_year, file.path(clean_path, paste0("occ_", y, ".parquet")))
  
  rm(occ_year, year_results); gc()
}
