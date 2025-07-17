# CDC Scenario Modeling COVID-19 Vaccine Coverage Data
# National Immunization Survey-Adult COVID Module (NIS-ACM)
# from here: 
# https://github.com/midas-network/covid19-scenario-modeling-hub/tree/main/auxiliary-data/vaccination-coverage/


# Only county-level data I've found in literature is from this disparities analysis tool
# https://data.cms.gov/tools/mapping-medicare-disparities-by-population
# Paper: https://www.ajpmonline.org/article/S0749-3797(21)00411-6/fulltext
# This is only for 65+ medicare vaccinated patients

#///////////////////////
#### LOAD LIBRARIES ####
#///////////////////////
library(tidyverse)
library(plotly)
library(tidycensus)

#/////////////////
#### ACS DATA ####
#/////////////////
source("../private_input_data/api_keys.R")
vars_acs5_2023 <- load_variables(2023, "acs5", cache = TRUE)

county_pops = 
  tidycensus::get_acs(geography = "county", state="TX", variables="B01001_001",
                      year = 2023, geometry=T)

#////////////////////////
#### VAX INTENT DATA ####
#////////////////////////
adult_vi = 
  read_csv("../../new_data_to_look_at/National_Immunization_Survey_Adult_COVID_Module__NIS-ACM___RespVaxView__Data___Centers_for_Disease_Control_and_Prevention__cdc.gov__20250707.csv")
ped_vi = 
  read.csv("../../new_data_to_look_at/Weekly_Parental_Intent_for_Vaccination_and_Cumulative_Percentage_of_Children_6_Months_-17_Years_Who_are_Up_to_date_with_the_COVID-19_Vaccines_by_Season__United_States_20250707.csv")

# Only sub state geographies for TX are
# "TX-Rest of State",   "TX-City of Houston", "TX-Bexar County" = San Antonio   
sub_adult_vi = adult_vi %>%
  filter(`Geography Type` == "Substate") %>%
  filter(Geography == "TX-Rest of State") %>%
  #filter(str_detect(Geography, "County")) %>%
  filter(`Group Name` == "Social Vulnerability Index (SVI) of county of residence") %>%
  filter(str_starts(`Indicator Name`, "Vaccination coverage")) %>%
  filter(`Suppression Flag` == 0) %>%
  separate(`95% CI (%)`, into = c("ci_lower", "ci_upper"), sep = " - ", convert = TRUE) %>%
  mutate(
    Time_Year = paste0(Year, " ", `Time Period`),
    # Extract first month and day
    first_month_day = str_extract(`Time Period`, "^[A-Za-z]+ \\d+"),
    # Combine with Year to form a date string
    first_date_str = paste0(first_month_day, " ", Year),
    # Parse to Date type (ensure locale is English if needed)
    first_date = mdy(first_date_str),
    # Social Vulnerability Index (SVI) of county of residence
    `Group Category` = factor(`Group Category`, levels = c("Low SVI", "Moderate SVI", "High SVI")),
    `Indicator Category` = 
      factor(`Indicator Category`, 
             levels=c("Received a 2024-2025 COVID-19 dose",
                      "Definitely will get the 2024-2025 COVID-19 vaccine",
                      "Probably will get the 2024-2025 COVID-19 vaccine or are unsure",
                      "Probably or definitely will not get the 2024-2025 COVID-19 vaccine"
                      )))
date_breaks = sort(unique(sub_adult_vi$first_date))
# Get second to last and last
#last_date = date_breaks[length(date_breaks)]
#second_last = date_breaks[length(date_breaks)-1]

vax_intent_adult = 
  ggplot(sub_adult_vi, 
         aes(x = first_date, y = `Estimate (%)`, 
             fill = `Group Category`, color = `Group Category`,
             text = paste0(
               "Estimate: ", `Estimate (%)`, "%\n",
               "95% CI: ", ci_lower, " - ", ci_upper, "%\n",
               "Sample Size: ", `Sample Size`
             ))) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                alpha = 1, position = position_dodge(width = 15)) +
  geom_point(size = 2, alpha = 1, position = position_dodge(width = 15)) +
  facet_wrap( ~ `Indicator Category`, ncol=1, scales = "free_y") + # 
  labs(x = "Date", y = "Estimate (%)") +
  scale_x_date(
    breaks = date_breaks,
    date_labels = "%b %d %Y"
  ) +
  # actually looks bad with the facet header and jittered points
  # ggbreak::scale_x_break(c(second_last, last_date))+
  scale_color_manual(values = c("#008080", "#EB7F86", "#B95E9A"))+
  scale_fill_manual( values = c("#008080", "#EB7F86", "#B95E9A"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=25, hjust=1),
        legend.position = "bottom")
 
ggplotly(vax_intent_adult, tooltip = "text") %>%
  layout(legend = list(orientation = "h", x = 0.5,           
                       xanchor = "center", y = -0.2))       

vax_intent_adult


#//////////////////////////
#### VAX COVERAGE DATA ####
#//////////////////////////
adult_flu_cov = read_csv("../../new_data_to_look_at/Weekly_Influenza_Vaccination_Coverage_and_Intent_for_Vaccination__Overall__by_Selected_Demographics_and_Jurisdiction__Among_Adults_18_Years_and_Older_20250707.csv") %>%
  rename(Estimate = Estimates)
ped_flu_cov   = read_csv("../../new_data_to_look_at/Weekly_Cumulative_Influenza_Vaccination_Coverage_by_Flu_Season__Selected_Demographics__and_Race_and_Ethnicity_Among_Children_6_Months-17_Years__United_States_20250707.csv")

# Need to standardize col names because some are title case and some are lower
names(adult_flu_cov) <- toupper(names(adult_flu_cov))
names(ped_flu_cov) <- toupper(names(ped_flu_cov))

# Bind rows
all_flu_cov_df = adult_flu_cov %>%
  bind_rows(ped_flu_cov) %>%
  separate(WEEK_ENDING, into=c("WEEK_ENDING_DATE", NA, NA), sep=" ") %>%
  mutate(WEEK_ENDING_DATE = as.Date(WEEK_ENDING_DATE, format = "%m/%d/%Y")) %>%
  rowwise() %>%
  mutate(UPPER_ESTIMATE = ESTIMATE + CI_HALF_WIDTH_95PCT,
         LOWER_ESTIMATE = ESTIMATE - CI_HALF_WIDTH_95PCT ) %>%
  ungroup()


age_groups = c("6 months-4 years", "5-17 years", "18-29 years", "30-39 years", "40-49 years", "50-64 years", "65+ years")
age_strat_all_flu_cov_df = all_flu_cov_df %>%
  filter(`GEOGRAPHIC LEVEL`=="National" & `DEMOGRAPHIC LEVEL`=="Age") %>%
  filter(INDICATOR_LABEL=="Up-to-date") %>%
  filter(INFLUENZA_SEASON %in% c("2023-2024", "2024-2025")) %>%
  filter(`DEMOGRAPHIC NAME` %in% age_groups) %>%
  #filter(INDICATOR_CATEGORY_LABEL == "Received a vaccination") %>%
  mutate(`DEMOGRAPHIC NAME` = factor(`DEMOGRAPHIC NAME`, levels=age_groups))



# "#1D271CFF" "#274637FF" "#2C715FFF" "#44A57CFF" "#819A7AFF" "#58A449FF" "#CEC917FF" 
#color_pal = rev(ghibli::ghibli_palettes$MarnieMedium2[1:7])
color_pal = c("#8C2B0E", "#C5692D", "#FEB359", "#B47E83", "#68434E", "#435F90", "#132F5B")
#color_pal = NatParksPalettes::natparks.pals("DeathValley")[1:7]
#color_pal = c("#f3e79b","#fac484","#f8a07e","#eb7f86","#ce6693","#a059a0","#5c53a5") # Sunset


age_vax_cov = ggplot(age_strat_all_flu_cov_df, 
       aes(x=WEEK_ENDING_DATE, y=ESTIMATE, group=`DEMOGRAPHIC NAME`))+
  geom_ribbon(aes(ymin = LOWER_ESTIMATE, ymax = UPPER_ESTIMATE, fill=`DEMOGRAPHIC NAME`), alpha=0.5) +
  geom_line(aes(color=`DEMOGRAPHIC NAME`))+
  geom_point(aes(color=`DEMOGRAPHIC NAME`))+
  facet_wrap(~INFLUENZA_SEASON, nrow=1, scales="free_x")+
  scale_x_date(
    date_labels = "%m/%d/%y"
  ) +
  labs(x="Ending Week", y="Vaccine Coverage Estimate",
       color="Age Group", fill="Age Group")+
  scale_color_manual(values = color_pal) +
  scale_fill_manual(values = color_pal) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "right" #"bottom"
    )

age_vax_cov

plotly::ggplotly(age_vax_cov) %>%
  layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2))


tx_all_flu_cov_df = all_flu_cov_df %>%
  filter(`GEOGRAPHIC LEVEL`=="State" & `GEOGRAPHIC NAME`=="Texas") %>%
  filter(INDICATOR_LABEL=="Up-to-date") %>%
  filter(INFLUENZA_SEASON %in% c("2023-2024", "2024-2025"))
  
tx_vax_cov = ggplot(tx_all_flu_cov_df, 
                     aes(x=WEEK_ENDING_DATE, y=ESTIMATE, group=`DEMOGRAPHIC NAME`))+
  geom_ribbon(aes(ymin = LOWER_ESTIMATE, ymax = UPPER_ESTIMATE, fill=`DEMOGRAPHIC NAME`), alpha=0.5) +
  geom_line(aes(color=`DEMOGRAPHIC NAME`))+
  geom_point(aes(color=`DEMOGRAPHIC NAME`))+
  facet_wrap(~INFLUENZA_SEASON, nrow=1, scales="free_x")+
  scale_x_date(
    date_labels = "%m/%d/%y"
  ) +
  labs(x="Ending Week", y="Vaccine Coverage Estimate",
       color="Age Group", fill="Age Group")+
  scale_color_manual(values = color_pal) +
  scale_fill_manual(values = color_pal) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "right" #"bottom"
  )


plotly::ggplotly(tx_vax_cov) %>%
  layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2))




#/////////////////////////////////////////
#### CMS MAPPING MEDICARE DISPARITIES ####
#/////////////////////////////////////////

mmd_county = read_csv("../../new_data_to_look_at/MMD_DATA/flu_vax_by_county_mapping-medicare-disparities_data.csv")




















