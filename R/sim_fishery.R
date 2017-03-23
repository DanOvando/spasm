sim_fishery <- function(fish, fleet, manager,num_patches = 10,sim_years = 25,...){

pop <-  expand.grid(year = 1:sim_years, patch = 1:num_patches, age = 1:fish$max_age) %>%
  mutate(numbers = NA, biomass = NA, ssb = NA, numbers_caught = NA) %>%
  as_data_frame() %>%
  arrange(year,patch,age)

n0_at_age <- fish$r0/num_patches * exp(-fish$m * ( 0:(fish$max_age - 1)))

n0_at_age[fish$max_age] <- n0_at_age[fish$max_age]  / (1 - exp(-fish$m))

b0_at_age <- n0_at_age * fish$weight_at_age

ssb0_at_age <- n0_at_age * fish$ssb_at_age

pop$numbers[pop$year == 1] <- rep(n0_at_age, num_patches)

pop <- pop %>%
  left_join(data_frame(age = 1:fish$max_age, ssb_at_age = fish$ssb_at_age,
                      weight_at_age = fish$weight_at_age), by = 'age')

y <- 1

model_phase <- 'burn'

for (y in 1:sim_years){

  # Move adults

  pop[pop$year == y,] <- move_adults(pop %>% filter(year == y))

  # grow and die


  # growth_and_death <-  pop[pop$year == y,] %>%
  #   group_by(patch) %>%
  #   mutate(g_and_d = list(grow_and_die(numbers, fish, fleet))) %>%
  #   ungroup() %>% {
  #   purrr::transpose(.$g_and_d)
  #   } %>%
  #   map(unlist)

  pop[pop$year == (y + 1),'numbers'] <- pop[pop$year == y,] %>%
    group_by(patch) %>%
    mutate(numbers = grow_and_die(numbers, fish, fleet)$survivors) %>%
    ungroup() %>%
    {
      .$numbers
    }

  pop[pop$year == y,'numbers_caught'] <- pop[pop$year == y,] %>%
    group_by(patch) %>%
    mutate(numbers_caught = grow_and_die(numbers, fish, fleet)$caught) %>%
    ungroup() %>%
    {
      .$numbers_caught
    }

  pop <- pop %>%
    mutate(ssb = numbers * ssb_at_age,
           biomass = numbers * weight_at_age,
           biomass_caught = numbers_caught * weight_at_age)

  # spawn

  pop$numbers[pop$year == (y + 1) &
                pop$age == 1] <-
    calculate_recruits(pop = pop[pop$year == (y + 1), ],
                       fish = fish,
                       num_patches = num_patches,
                       phase = model_phase)

  if (y == 1){

    fish$ssb0 <- pop %>%
      filter(year == 1) %>%
      summarise(ssb = sum(ssb)) %>%
      ungroup() %>%  {
        (.$ssb)
      }

    model_phase <- 'recruit'

  }

  print(y)

}

# pop %>%
#   group_by(year) %>%
#   summarise(num = sum(numbers)) %>%
#   ggplot(aes(year,num)) +
#   geom_point()

# pop %>%
#   group_by(year) %>%
#   summarise(ssb = sum(ssb), recs = sum(numbers[age == 1])) %>%
#   ggplot(aes(ssb, recs)) +
#   geom_point()

return(pop)

}