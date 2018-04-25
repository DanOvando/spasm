library(tidyverse)
library(FishLife)
library(spasm)

fish <-
  create_fish(
    scientific_name = "Lutjanus campechanus",
    query_fishlife = T,
    mat_mode = "length",
    time_step = 1,
    sigma_r = 0.1,
    price = 5,
    price_cv = 0,
    price_ac = 0.25,
    steepness = 0.6,
    r0 = 1000,
    rec_ac = 0.25
  )


fleet <- create_fleet(
  fish = fish,
  cost_cv =  0.25,
  cost_ac = 0.25,
  q_cv = 0.25,
  q_ac = 0.25,
  fleet_model = "open-access",
  theta = 0.5,
  cost = 2,
  sigma_effort = 0,
  length_50_sel = 0.25 * fish$linf,
  initial_effort = 0.1,
  profit_lags =  4,
  beta = 1.3
)

sim <- spasm::sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0),
  num_patches = 1,
  sim_years = 100,
  burn_year = 50,
  time_step = fish$time_step,
  est_msy = T,
  tune_costs = T,
  b_v_bmsy_oa = 0.75
)

sim %>%
  group_by(year) %>%
  summarise(b = sum(biomass),
            effort = unique(effort)) %>%
  ggplot(aes(year, b)) +
  geom_point()