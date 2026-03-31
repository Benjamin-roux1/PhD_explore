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

# Calculate deviation per species to global average RS 
mammals_range %>%
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