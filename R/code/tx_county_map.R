
#### LOAD LIBRARIES ####
library(tidyverse)
library(sf)
library(jsonlite)
library(tidycensus)
library(igraph)

#### GET COUNTY GEO ####
# ACS 2019-2023 TX County population map
source("../private_input_data/api_keys.R")
county_pops = 
  tidycensus::get_acs(geography = "county", state="TX", variables="B01001_001",
                      year = 2023, geometry=T)

#### TRAVEL DATA ####
county_fips = sort(county_pops$GEOID)
length(county_fips)

# 254 x 254 travel matrix
work_travel_df = read_csv("../../data/texas/work_matrix_rel.csv", col_names = F)
names(work_travel_df)    = county_fips

# Work travel on scale 0 - 1 so 
min(work_travel_df) # 0
max(work_travel_df) # 0.99

# Get county centroid coordinates
centroid_df <- county_pops %>%
  st_centroid() %>%
  mutate(CENTROID_X = st_coordinates(geometry)[,1],
         CENTROID_Y = st_coordinates(geometry)[,2]) %>%
  select(GEOID, CENTROID_X, CENTROID_Y) %>%
  st_drop_geometry()

# Get origin - dest pairs for county work travel
work_long <- work_travel_df %>%
  mutate(origin_FIPS = county_fips) %>%
  gather(key= "destination_FIPS", value = "travel_weight", -origin_FIPS) %>%
  filter(travel_weight > 0.0001) %>% # remove if weight between counties is 1%
  mutate(across(c(origin_FIPS, destination_FIPS), as.character)) %>%
  left_join(centroid_df, by=c("origin_FIPS"="GEOID")) %>%
  rename(origin_CENTROID_X = CENTROID_X, 
         origin_CENTROID_Y = CENTROID_Y) %>%
  left_join(centroid_df, by=c("destination_FIPS"="GEOID")) %>%
  rename(destination_CENTROID_X = CENTROID_X, 
         destination_CENTROID_Y = CENTROID_Y)

ggplot(work_long %>%
         filter(travel_weight <= 0.01), aes(x=travel_weight))+
  geom_histogram(binwidth=0.001, fill="white", color="black")+
  theme_bw()


#### OPEN JSON FILES####
# List all JSON files in folder
json_files <- list.files(
  path = "../../OUTPUT_small_deterministic_min1/",
  pattern = "OUTPUT_small_deterministic_\\d+\\.json$",
  full.names = TRUE
)

# Helper to extract day from filename
extract_day <- function(filename) {
  str_extract(filename, "\\d+(?=\\.json$)") %>% as.integer()
}

# Parse and flatten all files into a single data frame
deterministic_sim_df <- map_dfr(json_files, function(file) {
  day <- extract_day(file)
  sim_data <- fromJSON(file)$data
  
  data.frame(
    DAY_ID       = day,
    GEOID        = sim_data$fips_id,
    Susceptible  = sim_data$compartment_summary$S,
    Exposed      = sim_data$compartment_summary$E,
    Infect_Asymp = sim_data$compartment_summary$A,
    Treated      = sim_data$compartment_summary$T,
    Infect_Symp  = sim_data$compartment_summary$I,
    Recovered    = sim_data$compartment_summary$R,
    Deceased     = sim_data$compartment_summary$D
  ) 
})

#/////////////////////////////////////////////////////////////////////////
#### DIFF FIRST & LAST SIM ####
# Columns to select for pivot
infect_cols <- c("Susceptible", "Exposed", "Infect_Asymp", "Treated", "Infect_Symp", "Recovered", "Deceased")

# Get day min and max
day_min <- min(deterministic_sim_df$DAY_ID, na.rm = TRUE)
day_max <- max(deterministic_sim_df$DAY_ID, na.rm = TRUE)

# Filter and pivot for comparison
sim_diff <- deterministic_sim_df %>%
  filter(DAY_ID %in% c(day_min, day_max)) %>%
  select(GEOID, DAY_ID, all_of(infect_cols)) %>%
  pivot_wider(names_from = DAY_ID, values_from = all_of(infect_cols),
              names_sep = "_day_") %>%
  mutate(across(ends_with(paste0("_day_", day_max)),
                .fns = list(diff = ~ . - get(str_replace(cur_column(), paste0("_day_", day_max), paste0("_day_", day_min)))),
                .names = "{str_replace(.col, paste0('_day_', day_max), '')}_change"))


#/////////////////////////////////////////////////////////////////////////
#### JOIN GEOM ####
full_deterministic_sim_df = deterministic_sim_df %>%
  mutate(GEOID = str_pad(GEOID, width = 5, side = "left", pad = "0")) %>%
  rowwise() %>%
  mutate(infected = Exposed + Infect_Asymp + Infect_Symp) %>%
  ungroup() %>%
  left_join(county_pops, by="GEOID") %>%
  st_as_sf() # make sure geometry keeps its CRS, etc

#/////////////////////////////////////////////////////////////////////////
#### PLOT COUNTY ####
# Which counties started the sim with some exposed people
exposed_starting <- full_deterministic_sim_df %>%
  filter(DAY_ID == 0, Exposed > 0)

exposed_origins <- exposed_starting$GEOID

exposed_edges <- work_long %>%
  filter(origin_FIPS %in% exposed_origins) %>%
  filter(!(origin_FIPS == destination_FIPS))

county_infected_plot = 
  ggplot(full_deterministic_sim_df %>%
           filter(DAY_ID == day_max)) +
  geom_sf(aes(fill = Recovered, geometry=geometry), color = "white") +
  geom_sf( # draw black line around starting infectious counties
    data = exposed_starting,
    aes(geometry = geometry), fill=NA,
    size = 3, color = "black"
  ) +
  geom_segment( # Edge layer from exposed origins
    data = exposed_edges,
    aes(x = origin_CENTROID_X, y = origin_CENTROID_Y,
        xend = destination_CENTROID_X, yend = destination_CENTROID_Y,
        linewidth = travel_weight),
    color = "black", alpha = 0.1
  ) +
  scale_fill_gradientn(
    colors = rcartocolor::carto_pal(5, "Sunset"), name = "Recovered",
    labels = scales::comma
  ) +
  labs(x="", y="", linewidth="Travel Weight") + # title = "Population by County in CATRAC Region", 
  ggspatial::annotation_north_arrow(
    location = "bl",        # bottom left
    which_north = "true",   # "true" north
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.1, "in"),
    style = ggspatial::north_arrow_fancy_orienteering()
  ) +
  theme_minimal()

print(county_infected_plot)

ggsave(
  filename = paste0("../figures/recovered_county_map_deterministicSEATIRD_no-travel.png"),
  plot = county_infected_plot,
  bg = "white", 
  width = 6, height = 5, units = "in", dpi = 600
)










