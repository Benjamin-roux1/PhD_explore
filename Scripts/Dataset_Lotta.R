# Visualise the expert dataset
mountains <- st_read("GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp") 
expert <- readxl::read_excel("Lotta_Schultz/mountain_experts_validation.xlsx", sheet = 1)
mapview(mountains, zcol = "DBaseName", legend = FALSE)
expert %>%
  filter(validated == "yes") %>%
  anti_join(mountains_03, by = c("mountain_range" = "Level_03"))

# Keep only level 03
mountains_03 <- mountains %>%
  st_make_valid() %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) %>%
  mutate(Level_03 = coalesce(Level_03, Level_02)) %>% # some ranges doesn't have level03 so take level02 in such cases
  group_by(Level_01, Level_02, Level_03) %>%
  summarise(geometry = st_union(geometry))

ggplot() +
  geom_sf(data = mountains, color = "grey") +
  geom_sf(data = mountains_03, color = "red4", fill = "red") +
  theme_minimal() 

# Add all mountain ranges to the expert validation dataframe AND add column groups + validated for each new mountain range
all_combinations <- expand_grid(
  Level_03 = unique(mountains_03$Level_03),
  groups = c("mammals", "birds", "reptiles")
)

expert_mountain <- all_combinations %>%
  left_join(mountains_03, by = "Level_03") %>%        # add geometry
  left_join(expert, by = c("Level_03" = "mountain_range", "groups")) %>%  # add validated
  mutate(validated = ifelse(is.na(validated), "no", validated)) %>%
  select(Level_01, Level_02, Level_03, groups, validated, name, geometry) %>%
  st_as_sf()

expert_mountain %>%
  filter(groups == "reptiles" & validated == "yes") %>%
  ggplot() +
  geom_sf(data = mountains, color = "#d9d9d9", linewidth = 0.1) +
  geom_sf(fill = "#fd8d3c", color = "#e31a1c", linewidth = 0.1) +
  theme_minimal() +
  labs(
    title = "Expert-validated reptiles range"
  )

# Mapview of mountain ranges level 03
mapview(expert_mountain, zcol = "Level_03", legend = FALSE)

# Summary of dataset per continent
# continent centroids for bubble positions
continent_centroids <- tibble(
  Level_01 = unique(expert_mountain$Level_01),
  lon = c(20, 0, -30, 60, 70, -100, 135, -150, -60),
  lat = c(0, 80, 20, 50, -20, 45, -25, 0, -15)
)

# summarise validated per continent and group
summary_data <- expert_mountain %>%
  st_drop_geometry() %>%
  filter(validated == "yes") %>%
  group_by(Level_01, groups) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(continent_centroids, by = "Level_01") %>%
  mutate(
    # offset bubbles side by side per continent
    lon_offset = case_when(
      groups == "mammals" ~ lon - 8,
      groups == "birds"   ~ lon,
      groups == "reptiles"~ lon + 8
    )
  )

data("World")
# plot
ggplot() +
  geom_sf(data = World, fill = "grey90", color = "white", size = 0.1) +
  geom_point(data = summary_data, 
             aes(x = lon_offset, y = lat, size = n, color = groups),
             alpha = 0.7) +
  geom_text(data = summary_data,
            aes(x = lon_offset, y = lat, label = n),
            color = "white", size = 3, fontface = "bold") +
  scale_size_continuous(range = c(5, 20)) +
  scale_color_manual(values = c("mammals" = "#0570b0", "birds" = "#238b45", "reptiles" = "#e31a1c")) +
  labs(title = "Validated mountain ranges per group", color = "Group") +
  guides(size = "none") +
  theme_void() +
  theme(legend.position = "bottom")
