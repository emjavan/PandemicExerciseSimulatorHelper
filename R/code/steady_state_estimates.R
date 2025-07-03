#### LOAD LIBRARIES ####
library(tidyverse)
library(deSolve)

#### SEATIRD MODEL ####
SEATIRD <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N <- S + E + A + Treat + I + R + D
    
    dS <- -beta * (A + Treat + I) * S / N
    dE <- beta * (A + Treat + I) * S / N - tau * E
    dA <- tau * E - (kappa + gamma + nu) * A
    dT <- kappa * A - (chi + gamma + nu) * Treat
    dI <- chi * Treat - (gamma + nu) * I
    dR <- gamma * (A + Treat + I)
    dD <- nu * (A + Treat + I)
    
    list(c(dS, dE, dA, dT, dI, dR, dD))
  })
}

#### TIME ####
# Time vector — long enough to reach steady state
max_days = 500
times <- seq(0, max_days, by = 1)

#### NU ####
# Nu values by age and risk
nu_age_df = data.frame(
  AGE_GRP = c("0-4", "5-24", "25-49", "50-64", "65+"),
  nu_vect_low = c(0.0000223193, 0.0000409747, 0.0000837293, 0.0000618089, 0.0000089781)
) %>%
  mutate(nu_vect_high = nu_vect_low*9)

#### RISK ####
# Risk ratios uniformly applied to populations
risk_prop_df = data.frame(
  AGE_GRP = c("0-4", "5-24", "25-49", "50-64", "65+"),
  high_risk_ratios = c(0.0803804, 0.146264700734, 0.21454520795, 0.337052864254, 0.529643792782)
) %>%
  mutate(low_risk_ratios  = 1 - high_risk_ratios)

#### PARMAS ####
fixed_params = c(
  beta = 0.5,
  tau = 1/1.2,
  kappa = 1/1.9,
  chi = 1/1.0,
  gamma = 1/10
)

#### FULL DF ####
# Get populations supplied to the county sims by age and risk
# Travis,  Dallas,  Harris, El Paso,  Potter
county_fips_vect = c("453", "113", "201", "141", "375")
all_age_pops = read_csv("../../data/texas/county_age_matrix_small.csv")
full_age_risk_df = all_age_pops %>%
  filter(fips %in% county_fips_vect) %>%
  mutate(COUNTY_FIPS = paste0("48", fips)) %>%
  select(-fips) %>%
  gather(key="AGE_GRP", value="TOTAL_POP", -COUNTY_FIPS) %>%
  left_join(risk_prop_df, by="AGE_GRP") %>%
  rowwise() %>%
  mutate(LOW_RISK_POP = round(TOTAL_POP*low_risk_ratios, 0),
         HIGH_RISK_POP = TOTAL_POP - LOW_RISK_POP) %>%
  ungroup() %>%
  left_join(nu_age_df, by="AGE_GRP") %>%
  rename(
    RISK_POP_LOW = LOW_RISK_POP,
    RISK_POP_HIGH = HIGH_RISK_POP,
    RISK_RATIO_LOW = low_risk_ratios,
    RISK_RATIO_HIGH = high_risk_ratios,
    DEATH_RATE_NU_LOW = nu_vect_low,
    DEATH_RATE_NU_HIGH = nu_vect_high
  ) %>%
  pivot_longer(
    cols = c(RISK_POP_LOW, RISK_POP_HIGH,
             RISK_RATIO_LOW, RISK_RATIO_HIGH,
             DEATH_RATE_NU_LOW, DEATH_RATE_NU_HIGH),
    names_to = c(".value", "RISK_GRP"),
    names_pattern = "(.*)_(LOW|HIGH)"
  )  %>%
  dplyr::select(COUNTY_FIPS, AGE_GRP, RISK_GRP, everything())

#### LOOP OVER COUNTY-AGE-RISK ####
# Loop over age and risk to get all the sim runs to n days
inital_exposed_per_group = 10
full_df = data.frame() # initialize df
for (i in 1:nrow(full_age_risk_df)) {
  # Model initial conditions
  if(full_age_risk_df$AGE_GRP[i] == "25-49" & full_age_risk_df$RISK_GRP[i] == "LOW"){
    init <- c( # 48113, Dallas example in small stochastic/deterministic
      S = full_age_risk_df$RISK_POP[i] - inital_exposed_per_group,
      E = inital_exposed_per_group, # I'm assuming risk gets 5 of the 10 exposed??
      A = 0, Treat = 0, I = 0, R = 0, D = 0)
  }else{
    init <- c( # 48113, Dallas example in small stochastic/deterministic
      S = full_age_risk_df$RISK_POP[i],
      E = 0, A = 0, Treat = 0, I = 0, R = 0, D = 0 )
  } # end initial conditions by age group
  
  # Parameters
  params <- c( fixed_params,
    nu = full_age_risk_df$DEATH_RATE_NU[i]
  )
  
  # Solve the ODEs
  out <- ode(y = init, times = times, func = SEATIRD, parms = params,
             method = "euler", atol = 1e-8, rtol = 1e-8) # "lsoda"
  out_df <- as.data.frame(out) %>%
    mutate(
      COUNTY_FIPS   = full_age_risk_df$COUNTY_FIPS[i], # 5 counties
      AGE_GRP       = full_age_risk_df$AGE_GRP[i],     # 5 age groups
      RISK_GRP      = full_age_risk_df$RISK_GRP[i],    # 2 risk groups
      START_POP     = init["S"],
      START_EXP     = init["E"],
      DEATH_RATE_NU = full_age_risk_df$DEATH_RATE_NU[i])
  
  full_df = full_df %>%
    bind_rows(out_df)
} # end loop over age x risk group pop and nu by county

#### ODE SOLVE FINAL SOLUTION ####
# get rest of susceptibles for state of TX never exposed on our model
#  bc travel and age contact turned off
rest_of_tx_pop = all_age_pops %>%
  filter(!(fips %in% county_fips_vect)) %>%
  # mutate(fips = str_pad(fips, 3, "left", 0),
  #        COUNTY_FIPS = paste0("48", fips)) %>%
  select(-fips)
rest_of_tx_pop_total = sum(rest_of_tx_pop)

# Final solution
ode_solve_final = full_df %>%
  filter(time==max_days) %>%
  #group_by(COUNTY_FIPS, START_POP, START_EXP) %>%
  summarise(across(S:D, sum)) %>%
  mutate(S = S + rest_of_tx_pop_total)

#### R0 ESTIMATE PER GRP ####
# Expand by risk group
nu_df_long <- nu_age_df %>%
  tidyr::pivot_longer(cols = c(nu_vect_low, nu_vect_high),
                      names_to = "RISK_GRP",
                      values_to = "nu") %>%
  mutate(RISK_GRP = ifelse(RISK_GRP == "nu_vect_low", "LOW", "HIGH"))

# Compute R0 by group
r0_by_group <- nu_df_long %>%
  rowwise() %>%
  mutate(
    R0 = fixed_params["beta"] * (
      1 / (fixed_params["kappa"] + fixed_params["gamma"] + nu) +
        1 / (fixed_params["chi"] + fixed_params["gamma"] + nu) +
        1 / (fixed_params["gamma"] + nu)
    )) %>%
  ungroup()

print(r0_by_group)



















