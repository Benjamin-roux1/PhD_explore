# -----------
# library
library(terra); library(arrow); library(dplyr); library(ggplot2)
library(tidyr)
lapply(list.files("Functions/", pattern = "\\.R$", full.names = TRUE), source)

# -------
# Focus here on MAMMALS
# ---

# Import and filter Lotta's checklist
checklist <- read.csv("Lotta_files/vertebrate_data/vertebrate_data_Benjamin.csv")
mammals_GMBA <- checklist %>%
  filter(group == "mammals" & expert_validated == "yes") %>%
  distinct()

# Calculate range size
mammals_GMBA <- mammals_GMBA %>%
  mutate(GMBA_range = max_elevation - min_elevation)


# ------
# Import and process GBIF dataset
mammals_GBIF <- open_dataset("GBIF/data/Mammals_parquetclean")
# download DEM
dem <- terra::rast("elevation/demMountains_GLO90.tif")
# Collect Parquet dataset to R
mammals_GBIF <- mammals_GBIF %>%
  dplyr::select(species, genus, family, decimalLatitude, decimalLongitude, Level_02,
                Level_03) %>%
  collect()

# Select only occurrences coordinates
pts_mammals <- mammals_GBIF %>%
  dplyr::select(decimalLongitude, decimalLatitude)
# Add elevation
mammals_GBIF$elevation <- terra::extract(dem, pts_mammals)[,2]
# Calculate elevation max and min + range size for each species in each mountain range
mammals_GBIFrange <- mammals_GBIF %>%
  group_by(Level_02, Level_03, species) %>%
  summarise(
    NumberOcc = n(),
    GBIF_range = quantile(elevation, 0.95, na.rm = TRUE) - quantile(elevation, 0.05, na.rm = TRUE),
    GBIF_maxelev = quantile(elevation, 0.95, na.rm = TRUE),
    GBIF_minelev = quantile(elevation, 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(NumberOcc > 1)


# ------
# Compare Lotta's dataset with GBIF  
# First, compare the absolute number of species per mountain range/system
count_GMBA <- checklist %>%
  filter(group == "mammals") %>%
  distinct(Mountain_range, sciname) %>%
  count(Mountain_range, name = "n_checklist")

count_gbif <- mammals_GBIFrange %>%
  group_by(Level_02, Level_03) %>%
  summarise(
    n_gbif = n_distinct(species),
    n_occ = sum(NumberOcc, na.rm = TRUE), 
    mean_occ = mean(NumberOcc, na.rm = TRUE)
  )

comparison <- full_join(count_GMBA, count_gbif, 
                        by = c("Mountain_range" = "Level_03")) %>%
  mutate(difference = n_checklist - n_gbif) %>%
  filter(difference > 0)

comparison %>%
  ggplot(aes(x = Level_02, y = difference, fill = Level_02)) +
  geom_hline(yintercept =  0, col = "grey50", linetype = "dashed") +
  geom_col() +
  theme_perso() +
  theme(
    legend.position = "none"
      ) +
  labs(
    subtitle = "Difference of total number of species between GBMA and GBIF per mountain range",
    x = "Mountain ranges",
    y = "Difference of number of species (GMBA - GBIF)"
  )

# ------
# Compare range size
mammals_global <- semi_join(mammals_GMBA, mammals_GBIFrange, by = c("Mountain_range" = "Level_03", "sciname" = "species")) %>% # select only species found in Lotta's dataset
  left_join(mammals_GBIFrange, mammals_GMBA, by = c("Mountain_range" = "Level_03", "sciname" = "species"))

mammals_global <- mammals_global %>%
  rename(GMBA_minelev = min_elevation,
         GMBA_maxelev = max_elevation)

mammals_global <- mammals_global %>%
  tidyr::pivot_longer(cols = c(GMBA_range, GBIF_range,
                               GMBA_minelev, GMBA_maxelev,
                               GBIF_minelev, GBIF_maxelev),
                      names_to = c("source", ".value"),
                      names_pattern = "(GMBA|GBIF)_(.*)"
                      )

mammals_diff <- mammals_global %>%
  group_by(sciname, Mountain_range) %>%
  summarise(
    diff_range = range[source == "GMBA"] - range[source == "GBIF"],
    diff_maxelev = maxelev[source == "GMBA"] - maxelev[source == "GBIF"],
    diff_minelev = minelev[source == "GMBA"] - minelev[source == "GBIF"],
    NumberOcc = first(NumberOcc)
  )

mammals_diff %>%
  ggplot(aes(x = NumberOcc, y = diff_range)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey60") +
  geom_point() +
  geom_smooth(method = "gam") +
  scale_x_log10() +
  theme_perso()

mammals_diff %>%
  pivot_longer(cols = c(diff_maxelev, diff_minelev),
               names_to = "limit_elev",
               values_to = "diff_elev") %>%
  ggplot(aes(x = NumberOcc, y = diff_elev, colour = limit_elev, fill = limit_elev)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey60") +
  geom_point() +
  geom_smooth(method = "gam") +
  scale_x_log10() +
  theme_perso()




# ---------
# Boostrapping 
# Parameters
B <- 999
n_seq <- c(5, 10, 15, 20, 30, 40, 50, 60, 70, 85, 100, 120, 140, 160, 180, 200)
min_occ <- 200

# keep only species present in both datasets
mammals_GBIF <- mammals_GBIF %>%
  semi_join(mammals_global, by = c("species" = "sciname",
                                   "Level_03" = "Mountain_range"))

# Reference species from raw occurrences
reference_species <- mammals_GBIF %>%
  group_by(species, Level_02, Level_03) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= min_occ)

# True range from mammals_global (GMBA values)
true_ranges <- mammals_global %>%
  select(sciname, Mountain_range, GMBA_range, GMBA_maxelev, GMBA_minelev)

# Bootstrap function
bootstrap_stability <- function(data, true_range, true_maxelev, 
                                true_minelev, n_seq, B) {
  map_dfr(n_seq, function(n) {
    if (nrow(data) < n) return(NULL)
    
    replicates <- map_dfr(1:B, function(b) {
      sample_data <- slice_sample(data, n = n)
      tibble(
        ratio_range  = (quantile(sample_data$elevation, 0.95, na.rm = TRUE) -
                          quantile(sample_data$elevation, 0.05, na.rm = TRUE)) / 
          true_range,
        ratio_maxelev = (quantile(sample_data$elevation, 0.95, na.rm = TRUE) - 
          true_maxelev),
        ratio_minelev = (quantile(sample_data$elevation, 0.05, na.rm = TRUE) - 
          true_minelev)
      )
    })
    
    replicates %>%
      summarise(
        n                  = n,
        cv_range           = sd(ratio_range)   / mean(ratio_range),
        cv_maxelev         = sd(ratio_maxelev) / mean(ratio_maxelev),
        cv_minelev         = sd(ratio_minelev) / mean(ratio_minelev),
        mean_ratio_range   = mean(ratio_range,   na.rm = TRUE),
        mean_ratio_maxelev = mean(ratio_maxelev, na.rm = TRUE),
        mean_ratio_minelev = mean(ratio_minelev, na.rm = TRUE)
      )
  })
}

# Run bootstrap
library(furrr)
plan(multisession, workers = 10)

bootstrap_results <- reference_species %>%
  left_join(true_ranges, 
            by = c("species" = "sciname", 
                   "Level_03" = "Mountain_range")) %>%
  filter(!is.na(GMBA_range)) %>%
  future_pmap_dfr(function(species, Level_02, Level_03, n,
                           GMBA_range, GMBA_maxelev, GMBA_minelev) {
    data <- mammals_GBIF %>%
      filter(species == !!species, Level_03 == !!Level_03)
    
    bootstrap_stability(data, GMBA_range, GMBA_maxelev, 
                        GMBA_minelev, n_seq, B) %>%
      mutate(species  = species,
             Level_02 = Level_02,
             Level_03 = Level_03)
  }, .options = furrr_options(seed = 42))


bootstrap_results %>%
  group_by(n, Level_02) %>%
  summarise(
    ratio_range = mean(mean_ratio_range, na.rm = TRUE),
    ratio_maxelev = mean(mean_ratio_maxelev, na.rm = TRUE),
    ratio_minelev = mean(mean_ratio_minelev[is.finite(mean_ratio_minelev)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("ratio"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = n, y = value, color = metric)) +
  geom_hline(yintercept = 1, linewidth = 0.1, color = "grey40", linetype = "dashed") +
  geom_point() +
  geom_line(size = 1) +
  facet_wrap(~Level_02) +
  theme_perso() +
  theme(legend.position = "bottom")


bootstrap_results %>%
  group_by(n, Level_02) %>%
  summarise(
    ratio_range = mean(mean_ratio_range, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = ratio_range)) +
  geom_hline(yintercept = 1, linewidth = 0.1, color = "grey40", linetype = "dashed") +
  geom_point() +
  geom_line(size = 1) +
  facet_wrap(~Level_02) +
  theme_perso() +
  theme(legend.position = "bottom")
