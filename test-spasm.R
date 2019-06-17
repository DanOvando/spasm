library(tidyverse)
library(FishLife)
library(spasm)
library(ggridges)
library(sampling)
library(gganimate)

fish <-
  create_fish(
    scientific_name = "Atractoscion nobilis",
    query_fishlife = T,
    mat_mode = "length",
    time_step = 1,
    sigma_r = 0,
    price = 10,
    price_cv = 0,
    price_ac = 0,
    price_slope = 0,
    steepness = 0.9,
    r0 = 100,
    rec_ac = 0,
    density_movement_modifier = 0.1,
    adult_movement = 20,
    larval_movement = 3,
    density_dependence_form = 3
  )




fleet <- create_fleet(
  fish = fish,
  cost_cv =  0,
  cost_ac = 0,
  cost_slope = 0,
  q_cv = 0,
  q_ac = .7,
  q_slope = 0,
  fleet_model = "open-access",
  target_catch = 200,
  sigma_effort = 0,
  length_50_sel = 0.1 * fish$linf,
  initial_effort = 1000,
  profit_lags =  1,
  beta = 1,
  max_cr_ratio = 0.8,
  max_perc_change_f = 2,
  effort_allocation = 'profit-gravity',
  b_ref_oa = .9
)



a <- Sys.time()
sim_noad <- spasm::sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0.5, year_mpa = 5),
  num_patches = 50,
  sim_years = 50,
  burn_years = 1,
  time_step = fish$time_step,
  est_msy = FALSE,
  random_mpas = TRUE,
  min_size = 0.05,
  mpa_habfactor = 1,
  sprinkler = FALSE,
  keep_burn = FALSE,
  tune_costs = FALSE
)
Sys.time() - a
plot_spasm(sim_noad)
b = plot_spasm(sim_noad)
a = plot_spasm(sim_noad)

sim_noad %>%
  group_by(year, patch) %>%
  summarise(
    te = sum(effort),
    profits = sum(profits),
    biomass = sum(biomass),
    tc = sum(biomass_caught)
  ) %>%
  ungroup() %>%
  mutate(ppue = profits / te) %>%
  gather(metric, value,-year, -patch) %>%
  ggplot(aes(year, value, color = factor(patch))) +
  geom_line(show.legend = F) +
  facet_wrap( ~ metric, scales = "free_y")

sim_noad %>%
  group_by(year) %>%
  summarise(
    te = sum(effort),
    profits = sum(profits),
    biomass = sum(biomass),
    tc = sum(biomass_caught)
  ) %>%
  ungroup() %>%
  mutate(ppue = profits / te) %>%
  gather(metric, value,-year) %>%
  ggplot(aes(year, value)) +
  geom_line(show.legend = F) +
  facet_wrap( ~ metric, scales = "free_y")



sim_noad %>%
  group_by(year, patch) %>%
  summarise(
    biomass = sum(biomass),
    mpa = unique(mpa),
    effort = sum(effort)
  ) %>%
  # complete(biomass, nesting(year,patch), fill = list(biomass = 0)) %>%
  ungroup() %>%
  # filter(year == max(year)) %>%
  ggplot(aes(x = patch, y = effort, fill = mpa)) +
  geom_col(color = "transparent") +
  # geom_area(alpha = 0.5, na.rm = TRUE) +
  transition_time(year) +
  ease_aes('linear') +
  labs(title = 'Year: {frame_time}')


sim_noad %>%
  group_by(patch) %>%
  summarise(m = unique(eventual_mpa)) %>%
  ggplot(aes(patch, m)) +
  geom_point()

sim_noad %>%
  group_by(year, patch) %>%
  summarise(m = sum(ssb)) %>%
  ggplot(aes(year, m, color = factor(patch))) +
  geom_line(show.legend = FALSE)

sim_noad %>%
  group_by(year, patch) %>%
  summarise(ssb = sum(ssb),
            recs = numbers[age == 0]) %>%
  group_by(patch) %>%
  mutate(lssb = lag(ssb)) %>%
  ggplot(aes(lssb, recs)) +
  geom_line(color = "red") +
  facet_wrap( ~ patch)

sim_noad %>%
  group_by(year, patch) %>%
  summarise(biomass = sum(biomass),
            mpa = unique(mpa)) %>%
  # complete(biomass, nesting(year,patch), fill = list(biomass = 0)) %>%
  ungroup() %>%
  # filter(year == max(year)) %>%
  ggplot(aes(x = patch, y = biomass, fill = mpa)) +
  geom_col(color = "transparent") +
  # geom_area(alpha = 0.5, na.rm = TRUE) +
  transition_time(year) +
  ease_aes('linear') +
  labs(title = 'Year: {frame_time}')


sim_noad %>%
  mutate(q = f / effort) %>%
  ggplot(aes(year, q)) +
  geom_point()


sim_noad %>%
  group_by(year) %>%
  summarise(te = sum(effort)) %>%
  ungroup() %>%
  ggplot(aes(year, te)) +
  geom_line()

sim_noad %>%
  ggplot(aes(age, year, height = numbers, group = year)) +
  geom_density_ridges(stat = "identity") +
  labs(x = "Length (cm)", title = "Proportional Length Distribution")




fish <-
  create_fish(
    scientific_name = "Lutjanus campechanus",
    query_fishlife = T,
    mat_mode = "length",
    time_step = 1,
    sigma_r = 0,
    price = 5,
    price_cv = 0,
    price_ac = 0,
    steepness = 0.6,
    r0 = 1000,
    rec_ac = 0,
    density_movement_modifier = 1
  )


fleet <- create_fleet(
  fish = fish,
  cost_cv =  0.25,
  cost_ac = 0.25,
  q_cv = 0,
  q_ac = 0,
  fleet_model = "constant-effort",
  theta = 0.5,
  cost = 2,
  sigma_effort = 0,
  length_50_sel = 0.25 * fish$linf,
  initial_effort = 1000,
  profit_lags =  4,
  beta = 1.3
)

sim_ad <- spasm::sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0.5),
  num_patches = 10,
  sim_years = 100,
  burn_year = 50,
  time_step = fish$time_step,
  est_msy = F,
  tune_costs = F,
  b_v_bmsy_oa = 0.75,
  random_mpas = F
)

noad <- sim_noad %>%
  group_by(year, patch) %>%
  summarise(b = sum(biomass),
            effort = unique(effort),
            mpa = unique(mpa)) %>%
  mutate(adult_density_effect = "none")

ad <- sim_ad %>%
  group_by(year, patch) %>%
  summarise(b = sum(biomass),
            effort = unique(effort),
            mpa = unique(mpa)) %>%
  mutate(adult_density_effect = "full")


noad %>%
  bind_rows(ad) %>%
  ungroup() %>%
  filter(year > 140) %>%
  ggplot(aes(
    patch,
    year,
    height = b,
    group = interaction(year, adult_density_effect) ,
    fill = adult_density_effect
  )) +
  geom_density_ridges(stat = "identity", alpha = 0.5)