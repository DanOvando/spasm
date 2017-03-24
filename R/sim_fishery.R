#' \code{sim_fishery} simulates an age structured spatially explicit
#' model forward with fleets etc.
#'
#' @param fish
#' @param fleet
#' @param manager
#' @param num_patches
#' @param sim_years
#' @param ...
#'
#' @return a pop object with population and catch trajectories
#' @export
#'
#' @examples sim_fishery(fish = fish, fleet = fleet,...)
#'
sim_fishery <-
  function(fish,
           fleet,
           manager,
           num_patches = 10,
           sim_years = 25,
           burn_year = 10,
           ...) {
    pop <-
      expand.grid(
        year = 1:sim_years,
        patch = 1:num_patches,
        age = 1:fish$max_age
      ) %>%
      mutate(
        numbers = NA,
        biomass = NA,
        ssb = NA,
        numbers_caught = NA,
        mpa = F
      ) %>%
      as_data_frame() %>%
      arrange(year, patch, age)

    n0_at_age <-
      fish$r0 / num_patches * exp(-fish$m * (0:(fish$max_age - 1)))

    n0_at_age[fish$max_age] <-
      n0_at_age[fish$max_age]  / (1 - exp(-fish$m))

    b0_at_age <- n0_at_age * fish$weight_at_age

    ssb0_at_age <- n0_at_age * fish$ssb_at_age

    pop$numbers[pop$year == 1] <- rep(n0_at_age, num_patches)

    pop <- pop %>%
      left_join(
        data_frame(
          age = 1:fish$max_age,
          ssb_at_age = fish$ssb_at_age,
          weight_at_age = fish$weight_at_age
        ),
        by = 'age'
      )

    y <- 1

    model_phase <- 'burn'

    adult_move_grid <-
      expand.grid(source = 1:num_patches,
                  sink = 1:num_patches) %>%
      mutate(
        distance = source - sink,
        prob = 1 / ((2 * pi) ^ (1 / 2) * fish$adult_movement) * exp(-(distance) ^
                                                                      2 / (2 * fish$adult_movement ^ 2))
      ) %>%
      group_by(source) %>%
      mutate(prob_move = prob / sum(prob))

    adult_move_matrix <- adult_move_grid %>%
      ungroup() %>%
      select(source, sink, prob_move) %>%
      spread(sink, prob_move) %>%
      select(-source) %>%
      as.matrix()

    larval_move_grid <-
      expand.grid(source = 1:num_patches,
                  sink = 1:num_patches) %>%
      mutate(
        distance = source - sink,
        prob = 1 / ((2 * pi) ^ (1 / 2) * fish$adult_movement) * exp(-(distance) ^
                                                                      2 / (2 * fish$adult_movement ^ 2))
      ) %>%
      group_by(source) %>%
      mutate(prob_move = prob / sum(prob))

    larval_move_matrix <- larval_move_grid %>%
      ungroup() %>%
      select(source, sink, prob_move) %>%
      spread(sink, prob_move) %>%
      select(-source) %>%
      as.matrix()


    eventual_f <- fleet$eq_f

    fleet$eq_f <- 0

    for (y in 1:(sim_years - 1)) {
      # Move adults
      pop[pop$year == y &
            pop$age > 1,] <-
        move_fish(
          pop %>% filter(year == y, age > 1),
          fish = fish,
          num_patches = num_patches,
          move_matrix = adult_move_matrix
        )

      # change management

      if (y == manager$year_mpa) {
        prop_mpas <-  floor(num_patches * manager$mpa_size)

        mpa_locations <- sample(1:num_patches, prop_mpas)

        pop$mpa[pop$patch %in% mpa_locations] <-  T

        # reallocate fishing effort

        if (fleet$map_reaction == 'concentrate') {
          fleet$eq_f <- fleet$eq_f / (1 - length(mpa_locations) / num_patches)
        } else if (fleet$map_reaction == 'dilute') {
          fleet$eq_f <- fleet$eq_f * (1 - length(mpa_locations) / num_patches)

        }

      }
      # grow and die -----

      # growth_and_death <-  pop[pop$year == y,] %>%
      #   group_by(patch) %>%
      #   mutate(g_and_d = list(grow_and_die(numbers, fish, fleet))) %>%
      #   ungroup() %>% {
      #   purrr::transpose(.$g_and_d)
      #   } %>%
      #   map(unlist)

      pop[pop$year == (y + 1), 'numbers'] <-
        pop[pop$year == y,] %>%
        group_by(patch) %>%
        mutate(numbers = grow_and_die(numbers, mpa, fish, fleet)$survivors) %>%
        ungroup() %>%
        {
          .$numbers
        }

      pop[pop$year == y, 'numbers_caught'] <-
        pop[pop$year == y,] %>%
        group_by(patch) %>%
        mutate(numbers_caught = grow_and_die(numbers, mpa, fish, fleet)$caught) %>%
        ungroup() %>%
        {
          .$numbers_caught
        }

      pop <- pop %>%
        mutate(
          ssb = numbers * ssb_at_age,
          biomass = numbers * weight_at_age,
          biomass_caught = numbers_caught * weight_at_age
        )

      # spawn ----

      pop$numbers[pop$year == (y + 1) &
                    pop$age == 1] <-
        calculate_recruits(
          pop = pop[pop$year == (y + 1), ],
          fish = fish,
          num_patches = num_patches,
          phase = model_phase,
          move_matrix = larval_move_matrix
        )


      if (y == burn_year) {
        fish$ssb0 <- pop %>%
          filter(year == burn_year) %>%
          group_by(patch) %>%
          summarise(ssb = sum(ssb)) %>%
          ungroup() %>%  {
            (.$ssb)
          }

        model_phase <- 'recruit'

        fleet$eq_f <- eventual_f

      }


    }


    pop <- pop %>%
      filter(year > burn_year)

    return(pop)

  }