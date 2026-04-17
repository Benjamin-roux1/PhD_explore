library(terra); library(sf); library(leaflet); library(mapview)
library(ncdf4); library(arrow); library(dplyr); library(CopernicusDEM)
library(ggplot2); library(scattermore); library(tidyr)
library(sf); library(ggplot2); library(dplyr); library(tmap)
library(data.table)

plants <- open_dataset("GBIF/data/Tracheophyta_parquetclean")

# ---- Standardize species names
species_list <- plants %>%
  distinct(species) %>% 
  collect () %>%
  pull(species)
# use rgbif tool to find correct species names for each species in our dataset
matched <- rgbif::name_backbone_checklist(species_list)
# select the names that are not precise enough (genus, or less)
to_remove <- matched %>% 
  filter(matchType == "HIGHERRANK") %>% 
  pull(verbatim_name)
# Apply on the original Arrow dataset before collect()
plants_clean <- plants %>%
  filter(!species %in% to_remove)

# ---- Elevation dataset
dem <- terra::rast("elevation/demMountains_GLO90.tif")

# Calculate range size & number of occurrences
df_plants <- plants %>%
  dplyr::select(species, decimalLatitude, decimalLongitude, Level_03, DBaseName, Level_01,
                genus, family) %>%
  collect() %>%
  setDT()  # Use data.table instead of dplyr to handle very big dataset (+50M rows)

pts_plants <- df_plants[, .(decimalLongitude, decimalLatitude)]  # data.table syntax  
# spatvector better than 2 columns data frame because terra::extract expect this anyway
# BUT here it's too heavy

# Add elevation
# df_plants$elevation <- terra::extract(dem, pts_plants)[,2]
# OR do the following if df is too big (process by chunks)
chunk_size <- 1e6
df_plants[, elevation := NA_real_]  # data.table syntax + NA_real sets correct numerical type
for (i in seq(1, nrow(df_plants), by = chunk_size)) {
idx <- i:min(i + chunk_size - 1, nrow(df_plants))
df_plants[idx, elevation := terra::extract(dem, pts_plants[idx, ])[, 2]]
cat("Processed", max(idx), "/", nrow(df_plants), "\n")
}

# Calculate range size for each species in each mountain range
plants_range <- df_plants[
  !Level_01 %in% c("Arctic Ocean", "Atlantic Ocean", "Indian Ocean", "Pacific Ocean")
][,
  .(RangeSizeSpecies = quantile(elevation, 0.95, na.rm = TRUE) - quantile(elevation, 0.05, na.rm = TRUE),
    NumberOcc = .N),
  by = .(Level_03, species, Level_01, genus, family)
]

tab <- plants_range %>%
  mutate(occ_bin = cut(NumberOcc, breaks = c(0, 10, 100, 1000, Inf),
                       labels = c("<10", "10-100", "100-1000", ">1000"))) %>%
  count(Level_01, occ_bin)
# ggplot
ggplot(tab, aes(x = occ_bin, y = n, fill = occ_bin)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.4, size = 3, color = "grey30") +
  facet_wrap(~Level_01) +
  scale_fill_manual(values = c("#527E87FF", "#B88244FF", "#B8B69EFF", "#B48F2CFF")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey90", linewidth = 0.1)) +
  labs(
    title = "plants",
    x ="Number of occurrences",
    y = "Number of species"
  )
n_distinct(plants_range$species)

# range size ~ nb occurrences per continent
plants_range %>%
  filter(NumberOcc > 1) %>%
  slice_sample(n = 50000) %>%
  ggplot(aes(x = NumberOcc, y = RangeSizeSpecies)) +
  geom_point(alpha = 0.3, color = "blue4") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_x_log10() +
  theme_classic() +
  labs(
    title = "plants"
  )
# Mean range size per continent
plants_range %>%
  filter(NumberOcc > 1) %>%
  # slice_sample(n = 10000) %>%
  ggplot(aes(x = Level_01, y = RangeSizeSpecies, color = Level_01, fill = Level_01)) +
  geom_boxplot(alpha = 0.3) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  labs(
    title = "plants"
  )


# Calculate deviation to global average RS 
plants_range %>%
  filter(NumberOcc > 10) %>%
  mutate(global_mean = mean(RangeSizeSpecies, na.rm = TRUE),
         RS_dev = RangeSizeSpecies - global_mean) %>%
  ggplot(aes(x = Level_01, y = RS_dev, color = Level_01)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(
    y = "Deviation to mean global range size"
  )

# Calculate deviation per species to global average RS 
plants_range %>%
  filter(NumberOcc > 10) %>% #remove range size with too few occurrences
  group_by(species) %>% 
  filter(n_distinct(Level_03) > 1) %>% #remove species only found in ONE mountain range
  mutate(species_RS = mean(RangeSizeSpecies, na.rm = TRUE)) %>% #calculate global mean species range size
  ungroup() %>%
  mutate(SpeciesRS_dev = RangeSizeSpecies - species_RS) %>% #calculate deviation for each species in each mountain range
  ggplot(aes(x = Level_03, y = SpeciesRS_dev, color = Level_01)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~Level_01, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  labs(
    x = "Mountain range",
    y = "Range size deviation"
  )

# Link between deviation and number of occurrences
plants_range %>%
  filter(NumberOcc > 10) %>%
  group_by(species) %>%
  filter(n_distinct(Level_03) > 1) %>%
  mutate(species_RS = mean(RangeSizeSpecies, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(SpeciesRS_dev = RangeSizeSpecies - species_RS) %>%
  ggplot(aes(x = NumberOcc, y = SpeciesRS_dev, color = Level_01, fill = Level_01)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "loess", alpha = 0.1) + 
  scale_x_log10() +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey90", linewidth = 0.1)) +
  labs(
    x = "Number of Occurrences",
    y = "Range size deviation"
  )

# Link between RS and elevational distribution of the sampling effort
plants_all <- merge(df_plants, plants_range, by = c("species", "Level_03"), all.x = TRUE)  # data.table syntax of left_join

# Calculte elevational position of each occurrence relative to species range size
plants_elevpos <- plants_all[!is.na(elevation) & NumberOcc > 1]
plants_elevpos[, `:=` (
  elev_min = mean(elevation[elevation <= quantile(elevation, 0.05)], na.rm = TRUE),
  elev_max = mean(elevation[elevation >= quantile(elevation, 0.95)], na.rm = TRUE)
),
by = .(species, Level_03)
][, elevational_position := (elevation - elev_min) / (elev_max - elev_min)]

# Mean elevation position within the range function of the range size
# with every occurrences
plants_elevpos %>%
  filter(NumberOcc > 10) %>%
  ggplot(aes(x = RangeSizeSpecies, y = elevational_position)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE, interpolate = TRUE) +
  scale_fill_distiller(palette = "Greys", direction = 1, name = "Density") +
  geom_density_2d(color = "grey30", linewidth = 0.3, bins = 5) +
  geom_smooth(method = "lm", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 12) +
  labs(
    x = "Range size",
    y = "Relative elevational sampling position (all occurrences)"
  )

# with mean per species
plants_elevpos[, .(
  mean_position = mean(elevational_position, na.rm = TRUE),
  RangeSizeSpecies = first(RangeSizeSpecies),
  NumberOcc = first(NumberOcc)
), by = .(species, Level_03)][NumberOcc > 10] %>%
  ggplot(aes(x = RangeSizeSpecies, y = mean_position)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE, interpolate = TRUE) +
  scale_fill_distiller(palette = "Greys", direction = 1, name = "Density") +
  geom_density_2d(color = "grey30", linewidth = 0.3, bins = 10) +
  geom_smooth(method = "gam", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 12) +
  labs(
    x = "Range size",
    y = "Relative elevational sampling position (mean per species x MR)"
  )

# Mean elevation position within the range depending on the number of occurrences
plants_elevpos[, .(
  mean_position = mean(elevational_position, na.rm = TRUE),
  RangeSizeSpecies = first(RangeSizeSpecies),
  NumberOcc = first(NumberOcc)
), by = .(species, Level_03)] %>%
  ggplot(aes(x = NumberOcc, y = mean_position)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "gam", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  theme_classic() +
  scale_x_log10()

# Mean range size depending to the mean elevation of the species
plants_elevpos[!is.na(elevation) & NumberOcc > 10, .(
  mean_elev = mean(elevation, na.rm = TRUE),
  mean_RS = mean(RangeSizeSpecies, na.rm = TRUE)
), by = .(species, Level_03)] %>%
  ggplot(aes(x = mean_RS, y = mean_elev)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  theme_classic() 
