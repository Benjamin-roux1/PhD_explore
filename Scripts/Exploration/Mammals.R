library(terra); library(sf); library(leaflet); library(mapview)
library(ncdf4); library(arrow); library(dplyr); library(CopernicusDEM)
library(ggplot2); library(scattermore); library(tidyr)
library(sf); library(ggplot2); library(dplyr); library(tmap)

mammals <- open_dataset("GBIF/data/Mammals_parquetclean")

# ---- Standardize species names
species_list <- mammals %>%
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
mammals_clean <- mammals %>%
  filter(!species %in% to_remove)

# ---- Elevation dataset
dem <- terra::rast("elevation/demMountains_GLO90.tif")

# Calculate range size & number of occurrences
df_mammals <- mammals %>%
  dplyr::select(species, decimalLatitude, decimalLongitude, Level_03, DBaseName, Level_01,
                genus, family) %>%
  collect()

pts_mammals <- df_mammals %>%
  dplyr::select(decimalLongitude, decimalLatitude)

# Add elevation
df_mammals$elevation <- terra::extract(dem, pts_mammals)[,2]
# OR do the following if df is too big (process by chunks)
# chunk_size <- 500000
# df_Tracheophyta$elevation <- NA
# for (i in seq(1, nrow(df_Tracheophyta), by = chunk_size)) {
# idx <- i:min(i + chunk_size - 1, nrow(df_Tracheophyta))
# df_Tracheophyta$elevation[idx] <- terra::extract(dem, as.matrix(pts_Tracheophyta[idx, ]))[, 1]
# cat("Processed", max(idx), "/", nrow(df_Tracheophyta), "\n")
# }

mammals_range <- df_mammals %>%
  filter(!Level_01 %in% c("Arctic Ocean", "Atlantic Ocean", "Indian Ocean", "Pacific Ocean")) %>%
  group_by(Level_03, species, Level_01, genus, family) %>%
  summarise(RangeSizeSpecies = quantile(elevation, 0.95, na.rm = TRUE) - quantile(elevation, 0.05, na.rm = TRUE),
            NumberOcc = n(),
            .groups = "drop")

tab <- mammals_range %>%
  mutate(occ_bin = cut(NumberOcc, breaks = c(0, 10, 100, 1000, Inf),
                       labels = c("<10", "10-100", "100-1000", ">1000"))) %>%
  count(Level_01, occ_bin)
ggplot(tab, aes(x = occ_bin, y = n, fill = occ_bin)) +
  geom_col() +
  facet_wrap(~Level_01) +
  scale_fill_manual(values = c("#527E87FF", "#B88244FF", "#B8B69EFF", "#B48F2CFF")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey90", linewidth = 0.1)) +
  labs(
    title = "Mammals",
    x ="Number of occurrences",
    y = "Number of species"
  )

# range size ~ nb occurrences per continent
mammals_range %>%
  filter(NumberOcc > 1) %>%
  # slice_sample(n = 10000) %>%
  ggplot(aes(x = NumberOcc, y = RangeSizeSpecies)) +
  geom_point(alpha = 0.3, color = "blue4") +
  scale_x_log10() +
  facet_wrap(~Level_01) +
  theme_classic() +
  labs(
    title = "Mammals"
  )
# Mean range size per continent
mammals_range %>%
  filter(NumberOcc > 1) %>%
  # slice_sample(n = 10000) %>%
  ggplot(aes(x = Level_01, y = RangeSizeSpecies, color = Level_01, fill = Level_01)) +
  geom_boxplot(alpha = 0.3) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  labs(
    title = "Mammals"
  )


# Calculate deviation to global average RS 
mammals_range %>%
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

# Calculate deviation per species to global average RS (group by continent)
mammals_range %>%
  group_by(species) %>% 
  filter(NumberOcc > 10) %>%  # remove range size with too few occurrences
  filter(n_distinct(Level_03) > 1) %>%  # remove species only found in ONE mountain range
  mutate(species_RS = mean(RangeSizeSpecies, na.rm = TRUE)) %>%  # calculate global mean species range size
  ungroup() %>%
  mutate(SpeciesRS_dev = RangeSizeSpecies - species_RS) %>%  # calculate deviation for each species in each mountain range
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
# Calculate deviation per species to global average RS WITH LATITUDE
df_mammals %>%
  filter(!Level_01 %in% c("Arctic Ocean", "Atlantic Ocean", "Indian Ocean", "Pacific Ocean")) %>%
  group_by(Level_03, species, Level_01, genus, family) %>%
  summarise(RangeSizeSpecies = quantile(elevation, 0.95, na.rm = TRUE) - quantile(elevation, 0.05, na.rm = TRUE),
            NumberOcc = n(),
            mean_lat = mean(decimalLatitude, na.rm = TRUE),
            mean_lon = mean(decimalLongitude, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(species) %>% 
  filter(NumberOcc > 10) %>%  # remove range size with too few occurrences
  filter(n_distinct(Level_03) > 1) %>%  # remove species only found in ONE mountain range
  mutate(species_RS = mean(RangeSizeSpecies, na.rm = TRUE)) %>%  # calculate global mean species range size
  ungroup() %>%
  mutate(SpeciesRS_dev = RangeSizeSpecies - species_RS) %>%  # calculate deviation for each species in each mountain range
  ggplot(aes(x = mean_lat, y = SpeciesRS_dev)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  theme_classic() +
  labs(
    x = "Latitude",
    y = "Range size deviation"
  )

# Link between deviation and number of occurrences
mammals_range %>%
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
mammals_all <- left_join(df_mammals, mammals_range)

mammals_elevpos <- mammals_all %>%
  filter(!is.na(elevation)) %>%
  filter(NumberOcc > 1) %>%
  group_by(species, Level_03) %>%
  mutate(
    elev_min = min(elevation, na.rm = TRUE),
    elev_max = max(elevation, na.rm = TRUE),
    elevational_position = (elevation - elev_min) / (elev_max - elev_min)
    ) %>%
  ungroup()

# Mean elevation position within the range function of the range size
# with every occurrences
mammals_elevpos %>%
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
mammals_elevpos %>%
  filter(NumberOcc > 10) %>%
  group_by(species, Level_03) %>%
  summarise(
    mean_elevPos    = mean(elevational_position, na.rm = TRUE),
    RangeSizeSpecies = first(RangeSizeSpecies),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = RangeSizeSpecies, y = mean_elevPos)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE, interpolate = TRUE) +
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "Density") +
  geom_density_2d(color = "grey30", linewidth = 0.3, bins = 10) +
  geom_smooth(method = "lm", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 12) +
  labs(
    x = "Range size",
    y = "Relative elevational sampling position (mean per species x MR)"
  )

# Mean elevation position within the range function to the number of occurrences
mammals_elevpos %>%
  group_by(species, Level_03) %>%
  summarise(
    mean_position = mean(elevational_position, na.rm = TRUE),
    RangeSizeSpecies = first(RangeSizeSpecies),
    NumberOcc = first(NumberOcc)
  ) %>%
  ggplot(aes(x = NumberOcc, y = mean_position)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  theme_classic() +
  scale_x_log10()

# Mean range elevation ~ range size
mammals_all %>%
  filter(!is.na(elevation)) %>%
  filter(NumberOcc > 10) %>%
  group_by(species, Level_03) %>%
  summarise(
    mean_elev = mean(elevation, na.rm = TRUE),
    mean_RS = mean(RangeSizeSpecies, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = mean_RS, y = mean_elev)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  theme_classic() 

# Add elevation min/max of the mountain range
GMBA_mountains <- st_read("GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp")
GMBA_summary <- GMBA_mountains %>%
  st_drop_geometry() %>%
  group_by(Level_03) %>%
  summarise(
    Elev_Low  = min(Elev_Low),   # lowest point across sub-regions
    Elev_High = max(Elev_High)   # highest point across sub-regions
  )
mammals_all <- left_join(mammals_all, GMBA_summary, by = "Level_03")
# Species range size ~ highest elevation of the mountain range
mammals_all %>%
  filter(NumberOcc > 10) %>%
  mutate(Elev_range = Elev_High - Elev_Low) %>%
  group_by(species, Level_03) %>%
  summarise(
    mean_RS = mean(RangeSizeSpecies, na.rm = TRUE),
    Elev_range = first(Elev_range),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Elev_range, y = mean_RS)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  theme_classic() +
  labs(
    x = "Mountain elevational range",
    y = "Mean species range size"
  )

# Number of occurrences (sampling effort) ~ highest elevation MR
mammals_all %>%
  filter(NumberOcc > 10) %>%
  mutate(Elev_range = Elev_High - Elev_Low) %>%
  group_by(species, Level_03, Level_01) %>%
  summarise(
    mean_RS = mean(RangeSizeSpecies, na.rm = TRUE),
    Elev_range = first(Elev_range),
    NumberOcc = first(NumberOcc),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Elev_range, y = NumberOcc, color = Level_01, fill = Level_01)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "firebrick", fill = "firebrick", alpha = 0.15) +
  stat_ellipse(alpha = 0.1, geom = "polygon") +
  theme_classic() +
  scale_y_log10() +
  labs(
    x = "Mountain elevational range",
    y = "Number of occurrences"
  )

# Range size to mountain elevational range 
mammals_all %>%
  filter(NumberOcc > 10) %>%
  mutate(Elev_range = Elev_High - Elev_Low) %>%
  group_by(species, Level_03, Level_01) %>%
  summarise(
    Relative_RS = mean(RangeSizeSpecies, na.rm = TRUE)/mean(Elev_range, na.rm = TRUE),
    Elev_range = first(Elev_range),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Elev_range, y = Relative_RS, color = Level_01, fill = Level_01)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Level_01) + 
  theme_classic() +
  labs(
    x = "Mountain elevational range",
    y = "Relative range size"
  )


# ---------
