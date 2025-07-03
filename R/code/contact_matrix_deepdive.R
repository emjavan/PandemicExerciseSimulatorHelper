#///////////////////////////////////////////////////////////////////////////////
#' Compare contact matrices from Mistry 2021 and 
#'  what was provided by Lauren Ancel Meyers' (LAM's) group for model
#///////////////////////////////////////////////////////////////////////////////

library(tidyverse)
library(tidycensus)
source("../private_input_data/api_keys.R")

#### LAM GRP ####
# Contact matrix modified in transmission rate as the fraction of county population that is that age, risk, vax group
county_pop_by_age_matrix = read_csv("../../data/texas/county_age_matrix.csv") %>%
  mutate(fips = paste0("48", str_pad(fips, 3, "left", pad = "0"))) %>% # full county 5-digit fips
  rename(COUNTY_FIPS = fips)
age_grp_names = names(county_pop_by_age_matrix )[2:ncol(county_pop_by_age_matrix)]

county_pop_total = county_pop_by_age_matrix %>%
  rowwise() %>%
  mutate(TOTAL_COUNTY_POP = sum(c_across(`0-4`:`65+`), na.rm = TRUE)) %>%
  ungroup()

# Contacts are on a pretty huge scale
# Assume these have to be contacts with the other age group per day?? for the entire population and not a single person?
contact_mat_LAM = read_csv("../../data/texas/contact_matrix.csv", col_names = F)
names(contact_mat_LAM) = age_grp_names
min(contact_mat_LAM) # 4.022272
max(contact_mat_LAM) # 45.12285

contact_mat_LAM_long = contact_mat_LAM %>%
  mutate(AGE_GRP_FROM_ROW = age_grp_names) %>%
  gather(key="AGE_GRP_TO_COL", value="CONTACTS", -AGE_GRP_FROM_ROW)

# Risk ratio is proportion of the age group at increased risk of death by making the nu_vect *9
risk_ratios = read_csv("../../data/texas/high_risk_ratios.csv", col_names = F) %>%
  mutate(AGE_GRP = age_grp_names) %>%
  rename(HIGH_RISK_RATIO = X1)

county_pop_long = county_pop_total %>%
  gather(key="AGE_GRP", value="POPULATION", -COUNTY_FIPS, -TOTAL_COUNTY_POP) %>%
  left_join(risk_ratios, by="AGE_GRP") %>%
  rowwise() %>%
  mutate(HIGH_RISK_POP = round(HIGH_RISK_RATIO*POPULATION, 0),
         LOW_RISK_POP = POPULATION - HIGH_RISK_POP) %>%
  ungroup() %>%
  left_join(contact_mat_LAM_long, by=c("AGE_GRP"="AGE_GRP_FROM_ROW"), relationship = "many-to-many") %>%
  rename(AGE_GRP_FROM_ROW = AGE_GRP) %>%
  dplyr::select(COUNTY_FIPS, AGE_GRP_FROM_ROW, AGE_GRP_TO_COL, everything()) %>%
  rowwise() %>%
  mutate(
    # Proportion of the entire county population in this age and risk group
    HIGH_RISK_PROP = HIGH_RISK_POP/TOTAL_COUNTY_POP,
    LOW_RISK_PROP  = LOW_RISK_POP/TOTAL_COUNTY_POP,
    
    # The number of contacts given this age and risk groups' 
    HIGH_RISK_CONTACTS = HIGH_RISK_PROP*CONTACTS,
    LOW_RISK_CONTACTS = LOW_RISK_PROP*CONTACTS
  ) %>%
  ungroup()

#### Mistry 2021 ####
# Texas ALL = Work + School + Home + Community
# from https://github.com/epistorm/epydemix-data/blob/main/data/United_States_Texas/contact_matrices/mistry_2021/contacts_matrix_all.csv
# Check readme for matrix details https://github.com/mobs-lab/mixing-patterns
# M_ij measures the average number of contacts for an individual of age i with all of their contacts of age j
epydemix_tx_all_contact = read_csv("../public_input_data/Mistry_2021_contacts_matrix_all_TX_2021.csv", col_names = F)
#min(epydemix_tx_all_contact) # 0.007224034
#max(epydemix_tx_all_contact) # 3.982432
names(epydemix_tx_all_contact)  = seq(0, 84)
epydemix_tx_all_contact$AgeContactInitiator = seq(0, 84)

# Define age bins and labels
age_bins   = c(0, 5, 18, 50, 65, 86)
age_labels = c("0-4", "5-17", "18-49", "50-64", "65+")

epydemix_tx_all_contact_agegrp = epydemix_tx_all_contact %>%
  gather(key="AgeReceivingContact", value="MeanContacts", -AgeContactInitiator) %>%
  mutate(
    AgeContactInitiator = as.numeric(AgeContactInitiator),
    AgeReceivingContact = as.numeric(AgeReceivingContact),
    AgeGrp_ContactInitiator = cut(AgeContactInitiator, breaks = age_bins, labels = age_labels, right = FALSE),
    AgeGrp_ReceivingContact = cut(AgeReceivingContact, breaks = age_bins, labels = age_labels, right = FALSE)
  )

#### Mistry ALL HEAT MAP ####
contact_long <- epydemix_tx_all_contact %>%
  pivot_longer(
    cols = -AgeContactInitiator,
    names_to = "AgeContactee",
    values_to = "Contacts"
  ) %>%
  mutate(AgeContactee = as.integer(AgeContactee))  # convert from character to numeric

# Calculate label positions (midpoints)
label_positions <- tibble(
  age_group = age_labels,
  x = head(age_bins, -1),
  xend = tail(age_bins, -1)
  ) %>%
  mutate(midpoint = (x + xend - 1) / 2 - 2 )  # subtract 1 to keep integer midpoints

# Add this to your ggplot
mistry_contact_plt = ggplot(contact_long, aes(x = AgeContactee, y = AgeContactInitiator, fill = Contacts)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Contacts", option = "inferno", limits = c(0, 4)) +
  geom_vline(xintercept = age_bins - 0.5, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = age_bins - 0.5, color = "white", linewidth = 0.3) +
  # Add group labels on axes
  annotate("text", x = label_positions$midpoint, y = -3, label = label_positions$age_group, hjust = 0) +
  annotate("text", x = -3, y = label_positions$midpoint, label = label_positions$age_group, angle = 90, hjust = 0) +
  #coord_cartesian(xlim = c(0, 84), ylim = c(0, 84), clip = "off") +  # allow labels outside tiles
  labs(y="Age of Contact Initiator (i)", x="Age of Contactee (j)",
       title = "Per Capita Daily Contacts: TX (Mistry et al. 2021 - All)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(20, 20, 40, 80)  # bottom/left margin for labels
  )

ggsave(
  paste0("../figures/mistry_2021_percapita_contacts.png"),
  mistry_contact_plt,
  height=8, width=10, units="in", bg="white"
)

plotly::ggplotly(mistry_contact_plt)

#### Multiply by TX POP ####
acs_vars = tibble(acs_variable_code = sprintf("B01001_%0.3d", c(3:25, 27:49)), ## 3:25 males, 27:49 females
       age_grouping = rep(c(  '0-4',   '5-9', '10-14', '15-17', 
                            '18-19',    '20',    '21', '22-24',
                            '25-29', '30-34', '35-39', '40-44', 
                            '45-49', '50-54', '55-59', '60-61', 
                            '62-64', '65-66', '67-69', '70-74', 
                            '75-79', '80-84',  '85+'), 2))

tx_pop_by_age = get_acs(geography="state", state="TX", variables = acs_vars$acs_variable_code, 
                        geometry=FALSE, year = 2023) %>%
  left_join(acs_vars, by=c("variable"="acs_variable_code")) %>%
  group_by(GEOID, NAME, age_grouping) %>%
  summarise(estimate = sum(estimate), 
            moe = sum(moe))

# Final age groups we want
# "0-4"   "5-17"  "18-49" "50-64" "65+" 


# Disaggregate evenly into single years
# Assume uniform distribution of population for all ages within group
acs_expanded <- tx_pop_by_age %>%
  mutate(
    age_start = as.integer(str_extract(age_grouping, "^\\d+")),
    age_end = as.integer(str_extract(age_grouping, "\\d+$"))
  ) %>%
  replace_na(list(age_end = 85)) %>%
  rowwise() %>%
  mutate(
    age = list(seq(age_start, age_end))
  ) %>%
  unnest(age) %>%
  mutate(pop_single_year = estimate / (age_end - age_start + 1))

# Group 84 and 85+ to be 84+
final_pop_by_age_tx = acs_expanded %>%
  select(age, pop_single_year) %>%
  mutate(age = ifelse(age==85, 84, age)) %>%
  group_by(age) %>%
  summarise(pop_single_year = sum(pop_single_year)) %>%
  ungroup() %>%
  mutate(
    AgeGrp_ReceivingContact = cut(age, breaks = age_bins, labels = age_labels, right = FALSE),
    age = as.character(age)) %>%
  group_by(AgeGrp_ReceivingContact) %>%
  mutate(AgeGrp_TotalPop = sum(pop_single_year)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(PopWeightOfTotal = pop_single_year/AgeGrp_TotalPop) %>%
  ungroup()

epydemix_agegrp_popweight = epydemix_tx_all_contact_agegrp %>%
  select(-AgeGrp_ReceivingContact) %>%
  spread(AgeReceivingContact, MeanContacts) %>%
  group_by(AgeGrp_ContactInitiator) %>%
  summarise(across(`0`:`84`, sum), .groups = "drop") %>%
  gather(AgeReceivingContact, MeanContacts, -AgeGrp_ContactInitiator) %>%
  left_join(final_pop_by_age_tx, by=c("AgeReceivingContact"="age")) %>%
  rowwise() %>%
  mutate(PopWeightedContacts = MeanContacts*PopWeightOfTotal) %>%
  ungroup()

final_contact_mat_weighted = epydemix_agegrp_popweight %>%
  group_by(AgeGrp_ContactInitiator, AgeGrp_ReceivingContact) %>%
  summarise(PopWeightedContacts = sum(PopWeightedContacts)) %>%
  ungroup()

final_contact_mat_weighted_wide = final_contact_mat_weighted %>%
  spread(AgeGrp_ReceivingContact, PopWeightedContacts) %>%
  select(-AgeGrp_ContactInitiator)

# Copy put in data/texas/
write.table(final_contact_mat_weighted_wide,
            "../output_data_for_model/mistry_agegrp_contact_matrix.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)

#### Age Grp Heat Map ####
final_contact_mat_weighted <- final_contact_mat_weighted %>%
  mutate(
    AgeGrp_ReceivingContact = factor(AgeGrp_ReceivingContact, levels = age_labels),
    AgeGrp_ContactInitiator = factor(AgeGrp_ContactInitiator, levels = age_labels)
  )

grpd_contact_plt <- ggplot(final_contact_mat_weighted, 
                           aes(x = AgeGrp_ReceivingContact, y = AgeGrp_ContactInitiator, fill = PopWeightedContacts)) +
  geom_tile(color = "white") +  # adds white gridlines
  scale_fill_viridis_c(name = "Contacts", option = "inferno", limits = c(0, 12)) +
  labs(
    x = "Age of Contactee (j)", 
    y = "Age of Contact Initiator (i)",
    title = "Age Group Pop Weighted Daily Contacts: TX (Mistry et al. 2021 - All)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.margin = margin(20, 20, 40, 80)
  )

ggsave(
  paste0("../figures/mistry_2021_percapita_contacts_AgeGrpPopWeighted.png"),
  grpd_contact_plt,
  height=8, width=10, units="in", bg="white"
)


plotly::ggplotly(grpd_contact_plt)


















