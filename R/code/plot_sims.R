#///////////////////////////////////////////////////////////////////////////////
#' Compare stochastic and deterministic simulations
#' Set the folder paths appropriately for correct param sets
#///////////////////////////////////////////////////////////////////////////////


# Need to refactor this code to take in the option column names from figures
# i.e. did I vary vaccination, was is R0 or travel. 
# Need to generate figures in a robust way regardless of what params are varying

#///////////////////////
#### LOAD LIBRARIES ####
#///////////////////////
library(tidyverse)
library(jsonlite)
library(plotly)
library(tidycensus)
library(sf)
library(igraph)
library(scales)
library(gt)

#### GET COUNTY GEO ####
# ACS 2019-2023 TX County population map
source("../private_input_data/api_keys.R")
county_pops = 
  tidycensus::get_acs(geography = "county", state="TX", variables="B01001_001",
                      year = 2023, geometry=T)

#### TRAVEL DATA ####
county_fips = sort(county_pops$GEOID)

# 254 x 254 travel matrix of 0-1 proportion pop traveling
work_travel_df = read_csv("../../../PandemicExerciseSimulatorTACC/data/texas/work_matrix_rel.csv", 
                          col_names = F)
names(work_travel_df)    = county_fips

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

#///////////////////
#### HELPER FNS ####
#///////////////////
# Helper to extract day from filename
extract_day = function(filename) {
  str_extract(filename, "\\d+(?=\\.json$)") %>% as.integer()
} # end extract_day

extract_sim = function(filename) {
  str_extract(filename, "_sim\\d+") %>%
    str_remove("_sim") %>%
    as.integer()
} # end extract_sim

#///////////////////////////
#### LOOP OVER DATASETS ####
#///////////////////////////
# List all JSON files in folder
fips = c("48113", "48453", "48201", "48141", "48375")
# parent_dir_name = "simulations_with_NPIs/"
# dir_names_vect = c("NoTravel_WithinAgeGrpContactOnly_HighR0", "NoTravel_Mistry2021AllContact_HighR0", "Travel_Mistry2021AllContact_HighR0",
#                    "NoTravel_WithinAgeGrpContactOnly_LowR0",  "NoTravel_Mistry2021AllContact_LowR0",  "Travel_Mistry2021AllContact_LowR0")

parent_dir_name = "simulations_with_vax/"
dir_names_vect = c("Travel_Mistry2021AllContact_HighR0_Lag0", "Travel_Mistry2021AllContact_HighR0_Lag14")

fig_dir = paste0("../figures/", parent_dir_name)
if(!dir.exists(fig_dir)){
  dir.create(fig_dir)
}

summary_compare_df = data.frame()
runtime_df = data.frame()
for(i in 1:length(dir_names_vect)){
  #///////////////////
  #### JSON PATHS ####
  #///////////////////
  dir_name = dir_names_vect[i]
  dir_split = unlist(str_split(dir_name, pattern="_"))
  travel_param  = dir_split[1]
  contact_param = dir_split[2]
  r0_param      = dir_split[3]
  
  if(parent_dir_name == "simulations_with_vax/"){
    vax_param = dir_split[4]
  }
  
  stochastic_json_files <- list.files(
    path = paste0("../", parent_dir_name, dir_name, "/OUTPUT_small_stochastic_min1"),
    pattern = "\\.json$",
    full.names = TRUE,
    recursive = TRUE
  )
  deterministic_json_files <- list.files(
    path = paste0("../", parent_dir_name, dir_name, "/OUTPUT_small_deterministic_min1"),
    pattern = "OUTPUT_small_deterministic_\\d+\\.json$",
    full.names = TRUE
  )
  
  #////////////////////
  #### JSON to DFs ####
  #////////////////////
  # Parse and flatten all files into a single data frame
  stochastic_sim_df = map_dfr(stochastic_json_files, function(file) {
    sim = extract_sim(file)
    day = extract_day(file)
    sim_data = fromJSON(file)$data
    
    data.frame(
      SIM_ID       = sim,
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
  }) %>%
    mutate(GEOID = as.character(GEOID))
  
  # Parse and flatten all files into a single data frame
  deterministic_sim_df = map_dfr(deterministic_json_files, function(file) {
    day = extract_day(file)
    sim_data = fromJSON(file)$data
    
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
  }) %>%
    mutate(GEOID = as.character(GEOID))
  
  #//////////////
  #### PLOTS ####
  #//////////////
  # TX wide summary rather than county specific
  tx_lvl_epi_stochastic = stochastic_sim_df %>%
    #filter(GEOID %in% fips) %>%
    group_by(DAY_ID, SIM_ID) %>%
    summarise(across(Susceptible:Deceased, sum), 
              .groups = "drop")
  
  tx_lvl_epi_deterministic = deterministic_sim_df  %>%
    #filter(GEOID %in% fips) %>%
    group_by(DAY_ID) %>%
    summarise(across(Susceptible:Deceased, sum), 
              .groups = "drop")
  
  #////////////////////
  #### Susceptible ####
  #////////////////////
  sus_plt = ggplot()+
    geom_line(data = tx_lvl_epi_stochastic,
              aes(x=DAY_ID, y=Susceptible, group=SIM_ID), color="#FEB359", alpha=0.5, linewidth=0.7)+
    geom_line(data = tx_lvl_epi_deterministic,
              aes(x=DAY_ID, y=Susceptible), color="#132F5B", linewidth=0.7)+
    labs(x="Day")+
    scale_y_continuous(labels = scales::label_comma()) +
    theme_bw(base_size=25)
  
  ggsave(
    paste0(fig_dir, "susceptible_plt_", dir_name, ".png"),
    sus_plt,
    height=8, width=9, units="in", bg="white"
  )
  
  #////////////////////////////////
  #### Infectious Asymptomatic ####
  #////////////////////////////////
  inf_a_plt = ggplot()+
    geom_line(data = tx_lvl_epi_stochastic, 
              aes(x=DAY_ID, y=Infect_Asymp, group=SIM_ID), color="#C5692D" , alpha=0.5, linewidth=0.7)+
    geom_line(data = tx_lvl_epi_deterministic,
              aes(x=DAY_ID, y=Infect_Asymp), color="#132F5B", linewidth=0.7)+
    labs(x="Day", y="Infectious Asymp.")+
    scale_y_continuous(labels = scales::label_comma()) +
    theme_bw(base_size=25)
  
  ggsave(
    paste0(fig_dir, "infect_asymp_plt_", dir_name, ".png"),
    inf_a_plt,
    height=8, width=9, units="in", bg="white"
  )
  
  #//////////////////
  #### Recovered ####
  #//////////////////
  reco_plt = ggplot()+
    geom_line(data = tx_lvl_epi_stochastic,
              aes(x=DAY_ID , y=Recovered, group=SIM_ID), color="#B47E83", alpha=0.4, linewidth=0.7)+
    geom_line(data = tx_lvl_epi_deterministic,
              aes(x=DAY_ID , y=Recovered), color="#132F5B", linewidth=0.7)+
    labs(x="Day")+
    scale_y_continuous(labels = scales::label_comma()) +
    theme_bw(base_size=25)
  
  ggsave(
    paste0(fig_dir, "recovered_plt_", dir_name, ".png"),
    reco_plt,
    height=8, width=9, units="in", bg="white"
  )
  
  #/////////////////
  #### Deceased ####
  #/////////////////
  dead_plt = ggplot()+
    geom_line(data = tx_lvl_epi_stochastic,
              aes(x=DAY_ID , y=Deceased, group=SIM_ID), color="#8C2B0E", alpha=0.5, linewidth=0.7)+
    geom_line(data = tx_lvl_epi_deterministic,
              aes(x=DAY_ID , y=Deceased), color="#132F5B", linewidth=0.7)+
    labs(x="Day")+
    scale_y_continuous(labels = scales::label_comma()) +
    theme_bw(base_size=25)
  
  ggsave(
    paste0(fig_dir, "deceased_plt_", dir_name, ".png"),
    dead_plt,
    height=8, width=9, units="in", bg="white"
  )
  
  # interactive_dead_plt <- plotly::ggplotly(dead_plt)
  # tmp_file <- tempfile(fileext = ".html")
  # htmlwidgets::saveWidget(interactive_dead_plt, tmp_file, selfcontained = TRUE)
  # file.copy(tmp_file, paste0("../figures/deceased_plt_", dir_name, ".html"), overwrite=T)
  
  #/////////////////////////
  #### Combination Plot ####
  #/////////////////////////
  combo_plot = 
    cowplot::plot_grid(
      sus_plt   + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.margin = margin(10, 10, 10, 10)),  
      inf_a_plt + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.margin = margin(10, 10, 10, 10)), 
      reco_plt  + theme(plot.margin = margin(10, 10, 10, 10)), 
      dead_plt  + theme(plot.margin = margin(10, 10, 10, 10)), 
      nrow=2, align = "v", axis = "lr", rel_heights = c(0.9, 1))
  
  ggsave(
    paste0(fig_dir, "combo_plt_", dir_name, ".png"),
    combo_plot,
    height=8, width=12, units="in", bg="white"
  )
  
  #////////////////
  #### Exposed ####
  #////////////////
  #color_pal = ghibli::ghibli_palettes$MarnieMedium2[1:5]
  # color_pal is NatParksPalettes::natparks.pals("DeathValley")[1:7]
  color_pal = c("#8C2B0E", "#C5692D", "#FEB359", "#132F5B", "#435F90")
  filtered_det_df = deterministic_sim_df %>%
    filter(GEOID %in% fips)
  filtered_stoch_df = stochastic_sim_df %>%
    filter(GEOID %in% fips)
  
  max_days = max(c(filtered_stoch_df$DAY_ID, filtered_det_df$DAY_ID))
  common_x_breaks = pretty(c(0, max_days), n = 5)  # auto-tick spacing
  max_exposed = max(filtered_stoch_df$Exposed) - 0.4*max(filtered_stoch_df$Exposed) # move down 40% from the top
  mid_pt_x = 0.6*max_days
  stochastic_plot = ggplot()+
    geom_line(data = filtered_stoch_df,
              aes(x=DAY_ID, y=Exposed, group=interaction(GEOID, SIM_ID), color=GEOID), alpha=0.5, linewidth=0.7)+
    labs(x="Day", title = paste0("Stochastic (top) vs Deterministic (bottom)\n", dir_name))+
    annotate("label", x = mid_pt_x, y = max_exposed, label = "Initially Exposed\n10 Low risk, 25-49yr per County",
             size=5)+ # size in mm of height so about 14-15pt font
    scale_color_manual(values = color_pal)+
    scale_x_continuous(limits = c(0, max_days), breaks = common_x_breaks) +
    theme_bw(base_size = 14)+
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  
  det_plot = ggplot()+
    geom_line(data = filtered_det_df, 
              aes(x=DAY_ID , y=Exposed, group=GEOID, color=GEOID), linewidth=1)+
    labs(x="Simulation Day", color="County FIPS") +
    scale_color_manual(values = color_pal)+
    scale_x_continuous(limits = c(0, max_days), breaks = common_x_breaks) +
    theme_bw(base_size = 14)+
    theme(
      legend.position = "bottom"
    )
  
  final_exposed_compare_plt = 
    cowplot::plot_grid(stochastic_plot, det_plot, 
                       ncol=1, align = "v", axis = "lr")
  ggsave(
    paste0(fig_dir, "exposed_compare_plot_", dir_name, ".png"),
    final_exposed_compare_plt,
    height=8, width=9, units="in", bg="white"
  )
  
  #//////////////////////////////////
  #### TX MAP DETERMINISTIC PLOT ####
  #//////////////////////////////////
  # JOIN GEOM 
  full_deterministic_sim_df = deterministic_sim_df %>%
    mutate(GEOID = str_pad(GEOID, width = 5, side = "left", pad = "0")) %>%
    rowwise() %>%
    #mutate(infected = Exposed + Infect_Asymp + Infect_Symp) %>%
    ungroup() %>%
    left_join(county_pops, by="GEOID") %>%
    st_as_sf() # make sure geometry keeps its CRS, etc
  
  # Final day of simulation
  day_max = max(deterministic_sim_df$DAY_ID)
  
  # Create edges based on exposed start counties to their top travel
  exposed_starting <- full_deterministic_sim_df %>%
    filter(DAY_ID == 0, Exposed > 0)
  exposed_origins <- exposed_starting$GEOID
  
  county_reco_map = 
    ggplot(full_deterministic_sim_df %>%
             filter(DAY_ID == day_max)) +
    geom_sf(aes(fill = Recovered, geometry=geometry), color = "white") +
    geom_sf( # draw black line around starting infectious counties
      data = exposed_starting,
      aes(geometry = geometry), fill=NA, size = 3, color = "black"
    ) +
    scale_fill_gradientn( 
      colors = c("#F3E79B", "#FEB359", "#435F90", "#132F5B"), name = "Recovered",
      labels = scales::comma, limits = c(0, 40000), oob = scales::squish
    ) +
    labs(x="", y="", linewidth="Travel Weight") +
    ggspatial::annotation_north_arrow(
      location = "bl",        # bottom left
      which_north = "true",   # "true" north
      pad_x = unit(0.1, "in"),
      pad_y = unit(0.1, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    ) +
    theme_minimal()
 
  # Conditionally create exposed_edges only if it's a Travel scenario
  if (grepl("NoTravel", dir_name)) {
    exposed_edges <- NULL  # or tibble() if you need an empty df
  } else {
    exposed_edges <- work_long %>%
      filter(origin_FIPS %in% exposed_origins) %>%
      filter(!(origin_FIPS == destination_FIPS))
    
    county_reco_map = county_reco_map +
      geom_segment( # Edge layer from exposed origins
        data = exposed_edges,
        aes(x = origin_CENTROID_X, y = origin_CENTROID_Y,
            xend = destination_CENTROID_X, yend = destination_CENTROID_Y,
            linewidth = travel_weight),
        color = "black", alpha = 0.1
      )
  } # end if travel was used in model
  
  ggsave(
    filename = paste0(fig_dir, "recovered_county_map_", dir_name, ".png"),
    plot = county_reco_map,
    bg = "white", 
    width = 6, height = 5, units = "in", dpi = 600
  )
  
  #//////////////////////////////////
  #### DETERMINISTIC STEADYSTATE ####
  #//////////////////////////////////
  deterministic_final = deterministic_sim_df %>%
    filter(DAY_ID == max(DAY_ID))  %>%
    summarise(across(Susceptible:Deceased, sum), 
              MEDIAN_SIM_DAY = median(DAY_ID), # also the max/min/everything for the deterministic model
              .groups = "drop") %>%
    mutate(
      DATA_DIR = dir_name,
      SIM_TYPE = "DETERMINISTIC",
      VALUE = "converged_below_tolerance",
      TOTAL_SIMS = 1,
      TRAVEL_PARAM = travel_param,
      CONTACT_PARAM = contact_param,
      R0_PARAM = r0_param)
  
  summary_compare_df = summary_compare_df %>%
    bind_rows(deterministic_final)
  
  #///////////////////////////////
  #### STOCHASTIC STEADYSTATE ####
  #///////////////////////////////
  # start_pop_total_for_epi = stochastic_sim_df %>%
  #   filter(GEOID %in% fips) %>%
  #   filter(DAY_ID==0 & SIM_ID==1) %>%
  #   summarise(across(Susceptible:Deceased, sum))
  total_sims = length(unique(stochastic_sim_df$SIM_ID))
  end_pop_total_for_epi = stochastic_sim_df %>%
    group_by(SIM_ID) %>%
    filter(DAY_ID == max(DAY_ID)) %>%
    summarise(across(Susceptible:Deceased, sum), 
              MAX_SIM_DAY = max(DAY_ID),
              .groups = "drop") %>%
    mutate(
      DATA_DIR = dir_name,
      SIM_TYPE = "STOCHASTIC",
      TOTAL_SIMS = total_sims,
      TRAVEL_PARAM = travel_param,
      CONTACT_PARAM = contact_param,
      R0_PARAM = r0_param)
  
  mean_end_summary <- end_pop_total_for_epi %>%
    group_by(DATA_DIR, SIM_TYPE, TOTAL_SIMS, TRAVEL_PARAM, CONTACT_PARAM, R0_PARAM) %>%
    summarise(across(Susceptible:Deceased, mean), 
              MEDIAN_SIM_DAY = median(MAX_SIM_DAY),
              .groups = "drop") %>%
    mutate(VALUE = "mean")
  
  summary_compare_df = summary_compare_df %>%
    bind_rows(mean_end_summary)
  
  #/////////////////
  #### RUNTIMES ####
  #/////////////////
  single_runtime_df = 
    read_csv(paste0("../", parent_dir_name, dir_name, "/sim_runtime_stats.csv")) %>%
    rename(SIM_ID = sim,
           SIM_TYPE = model,
           RUNTIME_SEC = time_sec
           ) %>%
    mutate(SIM_TYPE = toupper(SIM_TYPE)) %>%
    left_join(deterministic_final %>%
                rename(MAX_SIM_DAY = MEDIAN_SIM_DAY) %>%
                bind_rows(end_pop_total_for_epi), by=c("SIM_TYPE", "SIM_ID")) %>%
    dplyr::select(-VALUE)
  
  runtime_df = bind_rows(runtime_df, single_runtime_df)

} # end loop over the different simulation sets

# Make the needed values rounded for display 
summary_compare_df = summary_compare_df %>%
  mutate(Susceptible = round(Susceptible, 2),
         Recovered   = round(Recovered, 2),
         Deceased    = round(Deceased, 2) )

#### SAVE CSV ####
write.csv(
  x=summary_compare_df,
  file="../sum_stats/sim_test_model_summary_compartments.csv",
  row.names = F
)

View(summary_compare_df)

#//////////////////////
#### RUN TIME PLOT ####
#//////////////////////
# make plot of sum(r+d) the run time to see how increase in people tracked slows down model
runtime_df2 = runtime_df %>%
  rowwise() %>%
  mutate(
         total_pop = sum(Susceptible, Recovered, Deceased),
         final_outbreak_size = sum(Recovered, Deceased),
         final_percent_infected = round((final_outbreak_size/total_pop)*100, 0)
         ) %>%
  ungroup() %>%
  mutate(R0_PARAM = ifelse(R0_PARAM=="HighR0", "High R0, ~5", "Low R0, ~1.5"),
         TRAVEL_PARAM = ifelse(TRAVEL_PARAM=="NoTravel", "No Travel", "Travel"),
         CONTACT_PARAM = ifelse(CONTACT_PARAM=="Mistry2021AllContact", "Full Contact", "Age Group Only Contact"), 
         SIM_TYPE = str_to_title(SIM_TYPE))

runtime_plot = 
  ggplot(runtime_df2, aes(x=MAX_SIM_DAY, y=RUNTIME_SEC))+
  geom_smooth(aes(group=DATA_DIR), method="lm", formula=y~x, se=FALSE, color="grey", size=1, show.legend=FALSE) +
  geom_point(aes(group=SIM_TYPE, color=SIM_TYPE, shape=SIM_TYPE, size=final_percent_infected), alpha=0.5)+
  facet_grid(TRAVEL_PARAM+CONTACT_PARAM~R0_PARAM, scales="free_y") +
  labs(x="Simulation Days", y="Run Time (sec)",
       shape="Model", color="Model", size="Percent Pop\nInfected") +
  theme_bw(base_size=18) +
  guides(
    shape = guide_legend(override.aes = list(size = 5)),
    color = guide_legend(override.aes = list(size = 5))
  )

ggsave(
  filename = "../figures/run_time_plot.png",
  plot = runtime_plot,
  bg = "white", 
  width = 12, height = 9, units = "in", dpi = 600
)


#////////////////////
#### DIFF MODELS ####
#////////////////////
group_vars = c("DATA_DIR", "SIM_TYPE", "TRAVEL_PARAM", "CONTACT_PARAM", "R0_PARAM")
median_runtime_df = runtime_df %>%
  group_by(!!!syms(group_vars)) %>%
  summarise(`Median Run Time (sec)` = median(RUNTIME_SEC), 
            .groups = "drop")

# Pull summary table out to match the table in power point
diff_models_df = summary_compare_df %>%
  dplyr::select(ends_with("_PARAM"), DATA_DIR, SIM_TYPE, Susceptible, Recovered, Deceased) %>%
  pivot_wider(
    names_from = SIM_TYPE,
    values_from = c(Susceptible, Recovered, Deceased)
  ) %>%
  mutate(
    diff_Susceptible = Susceptible_STOCHASTIC - Susceptible_DETERMINISTIC,
    diff_Recovered = Recovered_STOCHASTIC - Recovered_DETERMINISTIC,
    diff_Deceased = Deceased_STOCHASTIC - Deceased_DETERMINISTIC
  ) %>%
  # Create the diff rows
  mutate(SIM_TYPE = "DIFFERENCE") %>%
  dplyr::select(
    ends_with("_PARAM"), DATA_DIR, SIM_TYPE,
    Susceptible = diff_Susceptible,
    Recovered = diff_Recovered,
    Deceased = diff_Deceased
  ) %>%
  # Reshape the original deterministic and stochastic rows back to long
  bind_rows(
    summary_compare_df %>%
      select(ends_with("_PARAM"), DATA_DIR, SIM_TYPE, Susceptible, Recovered, Deceased, MEDIAN_SIM_DAY)
  ) %>%
  dplyr::select(DATA_DIR, R0_PARAM, CONTACT_PARAM, TRAVEL_PARAM, everything()) %>%
  left_join(median_runtime_df, by=group_vars) %>%
  mutate(
    CONTACT_PARAM = factor(CONTACT_PARAM, levels = c("WithinAgeGrpContactOnly", "Mistry2021AllContact")),
    SIM_TYPE = factor(SIM_TYPE, levels = c("STOCHASTIC", "DETERMINISTIC", "DIFFERENCE"))
  ) %>%
  group_by(DATA_DIR) %>%
  arrange(R0_PARAM, CONTACT_PARAM, TRAVEL_PARAM, SIM_TYPE) %>%
  mutate( # 2 different ways, first assumes correct ordering has already happened and is less desirable 
    MEDIAN_SIM_DAY = 
      ifelse(SIM_TYPE == "DIFFERENCE",
             lag(MEDIAN_SIM_DAY, 2) - lag(MEDIAN_SIM_DAY, 1),
             MEDIAN_SIM_DAY),
    `Median Run Time (sec)` = 
      ifelse(SIM_TYPE == "DIFFERENCE",
             `Median Run Time (sec)`[SIM_TYPE == "STOCHASTIC"] - `Median Run Time (sec)`[SIM_TYPE == "DETERMINISTIC"],
             `Median Run Time (sec)`)
  ) %>%
  ungroup()
  
  
write.csv(
  x=diff_models_df,
  file="../sum_stats/sim_test_model_summary_compartments_diff.csv",
  row.names = F
)

#/////////////////////
#### PRETTY TABLE ####
#/////////////////////
pretty_diff_models_df = diff_models_df %>%
  # Clean up R0_PARAM to show only High or Low
  mutate(R0 = ifelse(R0_PARAM == "HighR0", "High",
                     ifelse(R0_PARAM == "LowR0", "Low", R0_PARAM)),
         # Format numeric columns to comma-separated with 2 decimals
         Susceptible = scales::comma(Susceptible, accuracy = 0.01),
         Recovered = scales::comma(Recovered, accuracy = 0.01),
         Deceased = scales::comma(Deceased, accuracy = 0.01),
         CONTACT_PARAM = ifelse(CONTACT_PARAM=="Mistry2021AllContact", "Mistry et al.\n2021", "Within Age Grp"),
         TRAVEL_PARAM = ifelse(TRAVEL_PARAM=="NoTravel", "No", "Yes"),
         SIM_TYPE = str_to_title(SIM_TYPE)
  ) %>%
  # Rename columns for presentation clarity
  rename(
    `Contact Pattern` = CONTACT_PARAM,
    Travel = TRAVEL_PARAM,
    `Simulation Type` = SIM_TYPE,
    `Median Sim Days` = MEDIAN_SIM_DAY
  ) %>%
  # Select and order desired columns
  dplyr::select(R0, `Contact Pattern`, Travel, `Simulation Type`,
         Susceptible, Recovered, Deceased, `Median Run Time (sec)`, `Median Sim Days`)

pretty_table = pretty_diff_models_df %>%
  group_by(R0, `Contact Pattern`, Travel) %>%
  mutate(
    R0 = ifelse(row_number() == 1, R0, ""),
    `Contact Pattern` = ifelse(row_number() == 1, `Contact Pattern`, ""),
    Travel = ifelse(row_number() == 1, Travel, "")
  ) %>%
  ungroup() %>%
  gt() %>%
  # Bold header
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  # Bold DIFFERENCE rows
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = `Simulation Type` == "Difference"
    )) %>%
  # Add table options for spacing
  tab_options(
    table.font.size = 18,
    data_row.padding = px(4),
    table.width = px(1200)
  )

pretty_table

gtsave(pretty_table,
       "../sum_stats/diff_models_clean_table.png",
       vwidth = 2000,  # table width in pixels
       expand = 0)     # no extra padding













