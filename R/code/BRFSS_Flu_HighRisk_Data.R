# BRFSS DATA!
# https://www.cdc.gov/brfss/annual_data/annual_data.htm
# Want to understand fraction of high risk individuals who received a FLU vaccine in last year

# Load libraries
library(tidyverse)

# Read .XPT file using haven
brfss <- haven::read_xpt("../public_input_data/BRFSS/LLCP2023.XPT") # .XPT is the SAS transport file

colnames(brfss)[grep("AGE", colnames(brfss), ignore.case = T)]

# Filter for Texas respondents
tx <- brfss %>%
  filter(`_STATE` == 48) %>% 
  # _AGEG5YR = age reported in 5 year age categories
  # FLUSHOT7 = Adult flu shot/spray in last 12 months
  # BPHIGH6 = blood pressure high
  # there are more high risk columns but not pulling for now (would need to do Remy's PLACES combination??)
  dplyr::select(`_STATE`, `_AGEG5YR`, FLUSHOT7, BPHIGH6, `_LLCPWT`) %>%
  rename(AGE_5YR_CAT = `_AGEG5YR`) %>%
  mutate(
    flu_shot = case_when(
      FLUSHOT7 == 1 ~ 1, 
      FLUSHOT7 == 2 ~ 0,
      FLUSHOT7 %in% c(7,9) ~ NA_real_
    ),
    flu_shot = factor(flu_shot),
    high_risk = case_when(
      BPHIGH6 == 1 ~ 1, # Yes
      BPHIGH6 == 2 ~ 0, # No
      BPHIGH6 %in% c(7,9) ~ NA_real_ # Not sure or Refused
    ),
    high_risk = factor(high_risk),
    age_group = case_when(
      AGE_5YR_CAT == 1 ~ "18-24",
      AGE_5YR_CAT == 2 ~ "25-29",
      AGE_5YR_CAT == 3 ~ "30-34",
      AGE_5YR_CAT == 4 ~ "35-39",
      AGE_5YR_CAT == 5 ~ "40-44",
      AGE_5YR_CAT == 6 ~ "45-49",
      AGE_5YR_CAT == 7 ~ "50-54",
      AGE_5YR_CAT == 8 ~ "55-59",
      AGE_5YR_CAT == 9 ~ "60-64",
      AGE_5YR_CAT == 10 ~ "65-69",
      AGE_5YR_CAT == 11 ~ "70-74",
      AGE_5YR_CAT == 12 ~ "75-79",
      AGE_5YR_CAT == 13 ~ "80+",
      AGE_5YR_CAT == 14 ~ NA_character_,  # Don't know / Refused / Missing
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group)
  ) %>%
  drop_na() %>% # 3,706
  mutate(
    risk_grp = if_else(high_risk == 1, "High Risk", "Low Risk"),
    flu_vax = if_else(flu_shot == 1, "Vaccinated", "Unvaccinated"),
    risk_vax_grp = paste(risk_grp, flu_vax, sep = "_")
  )

# Proportion of high risk individuals who received a flu shot
prop_by_age <- tx %>%
  group_by(age_group, risk_vax_grp) %>%
  summarise(
    total_in_group = n()
  ) %>%
  group_by(age_group) %>%
  mutate(
    prop_in_group = total_in_group / sum(total_in_group)
  ) %>%
  ungroup() %>%
  arrange(age_group) %>%
  mutate( # gives some coersion warnings that don't matter
    age_midpt = case_when(
      grepl("\\+", age_group) ~ 82,  # adjust as appropriate
      TRUE ~ as.numeric(sub("-.*", "", age_group)) + 
        (as.numeric(sub(".*-", "", age_group)) - as.numeric(sub("-.*", "", age_group))) / 2
    )
  )

#### PLOT DATA ####
ggplot(prop_by_age, aes(x=age_midpt, y=prop_in_group, group=risk_vax_grp, color=risk_vax_grp)) +
  geom_line() +
  geom_point(aes(size=total_in_group)) +
  #geom_smooth(method="lm", se=FALSE) +
  scale_x_continuous(breaks = prop_by_age$age_midpt, labels = prop_by_age$age_group) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_color_manual(
    values = c(
      "High Risk_Vaccinated" = "#8C2B0E", 
      "High Risk_Unvaccinated" = "#C5692D",
      "Low Risk_Vaccinated" =   "#132F5B",
      "Low Risk_Unvaccinated" = "#435F90"
    )) +
  labs(x = "Age Group", y="Proportion in Group", color="Risk and Vax", size="Sample Size") + 
  #geom_hline(yintercept = 0.7, linetype = "dashed") +
  # annotate("text",
  #          x = min(prop_by_age$age_midpt),  # left-most position
  #          y = 0.7,
  #          label = "Healthy Humans Target",
  #          vjust = -0.5,
  #          hjust = 0,
  #          size = 4) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )
 
 












