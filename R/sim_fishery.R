#' \code{sim_fishery} simulates an age structured spatially explicit
#' model forward with fleets etc.
#'
#' @param fish
#' @param fleet
#' @param manager
#' @param num_patches
#' @param sim_years
#' @param ...
#' @param burn_year
#' @param crashed_pop
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
           sim_years = 1,
           burn_year = 10,
           crashed_pop = 1e-3,
           ...) {

    sim_years <- burn_year + sim_years

    if (fleet$fleet_model == 'supplied-catch'){

      sim_years <- sim_years + 1

    }


    pop <-
      expand.grid(
        year = 1:sim_years,
        patch = 1:num_patches,
        age = fish$min_age:fish$max_age
      ) %>%
      mutate(
        numbers = NA,
        biomass = NA,
        ssb = NA,
        numbers_caught = NA,
        profits = NA,
        effort = 0,
        f = 0,
        mpa = F,
        cost = fleet$cost
      ) %>%
      as_data_frame() %>%
      arrange(year, patch, age)

    if (fleet$cost_function == 'distance from port'){

      cost_frame <- data_frame(patch = 1:num_patches) %>%
        mutate(cost = fleet$cost * (1 + fleet$cost_slope * (patch - 1)))

      pop <- pop %>%
        select(-cost) %>%
        left_join(cost_frame, by = 'patch')

    }

    effort <- vector(mode = 'double', length = sim_years)

    f <- vector(mode = 'double', length = sim_years)

    q <- rep(fleet$q, sim_years)

    rec_devs <- rnorm(sim_years, mean = -(fish$sigma_r ^ 2) / 2, sd = fish$sigma_r)
    ## autocorrelated recruitment deviations
    for (t in 2:length(rec_devs)) {
      rec_devs[t] <-
        rec_devs[t - 1] * fish$rec_ac + sqrt(1 - fish$rec_ac ^ 2) * rec_devs[t]
    }

    mpa_locations <- -1

    n0_at_age <-
      fish$r0 / num_patches * exp(-fish$m * (fish$min_age:fish$max_age))

    n0_at_age[fish$max_age] <-
      n0_at_age[fish$max_age]  / (1 - exp(-fish$m))

    b0_at_age <- n0_at_age * fish$weight_at_age

    ssb0_at_age <- n0_at_age * fish$ssb_at_age

    pop$numbers[pop$year == 1] <- rep(n0_at_age, num_patches)

    pop <- pop %>%
      left_join(
        data_frame(
          age = fish$min_age:fish$max_age,
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
        prob = 1 / ((2 * pi) ^ (1 / 2) * fish$larval_movement) * exp(-(distance) ^
                                                                      2 / (2 * fish$larval_movement ^ 2))
      ) %>%
      group_by(source) %>%
      mutate(prob_move = prob / sum(prob))

    larval_move_matrix <- larval_move_grid %>%
      ungroup() %>%
      select(source, sink, prob_move) %>%
      spread(sink, prob_move) %>%
      select(-source) %>%
      as.matrix()


    # eventual_f <- fleet$eq_f
    #
    # fleet$eq_f <- 0
    for (y in 1:(sim_years - 1)) {
      # Move adults

      now_year <- pop$year == y

      pop[now_year &
            pop$age > fish$min_age,] <-
        move_fish(
          pop %>% filter(year == y, age > fish$min_age),
          fish = fish,
          num_patches = num_patches,
          move_matrix = adult_move_matrix
        )

      # change management

      if ((y - burn_year) == manager$year_mpa) {
        prop_mpas <-  floor(num_patches * manager$mpa_size)

        mpa_locations <- sample(1:num_patches, prop_mpas)

        pop$mpa[pop$patch %in% mpa_locations & pop$year >= y] <-  T

        # reallocate fishing effort
        # if (fleet$mpa_reaction == 'concentrate' &
        #     fleet$fleet_model == 'constant-effort') {
        #   # fleet$eq_f <- fleet$eq_f / (1 - length(mpa_locations) / num_patches)
        #   # effort[y] <-
        #   #   effort[y - 1] / (1 - length(mpa_locations) / num_patches)
        #
        # }

      }
      # fleet response

      if (fleet$fleet_model == 'constant-catch' & y > burn_year){

        effort_for_catch <- nlminb(1, catch_target, target_catch =
                 fleet$target_catch,
               pop = pop %>% filter(year == y),
               num_patches = num_patches,
               mpa = mpa,
               fleet = fleet,
               lower = 0,
               use = 'opt',
               fish = fish)

        effort[y] <- effort_for_catch$par
        # popcheck <- catch_target(total_effort = effort_for_catch$par,target_catch =
        #                fleet$target_catch,
        #              pop = pop %>% filter(year == y),
        #              num_patches = num_patches,
        #              mpa = mpa,
        #              fleet = fleet,
        #              use = 'check',
        #              fish = fish) %>% sum()
        #
        # if (popcheck < crashed_pop){
        #
        #   warning('constant catch killing population')
        #
        # }

      }

      if (fleet$fleet_model == 'supplied-catch' & y > burn_year){

        target_catch <- fleet$catches[y - burn_year]

        effort_for_catch <- nlminb(1, catch_target, target_catch =
                                     target_catch,
                                   pop = pop %>% filter(year == y),
                                   num_patches = num_patches,
                                   mpa = mpa,
                                   fleet = fleet,
                                   lower = 0,
                                   use = 'opt',
                                   fish = fish)

        effort[y] <- effort_for_catch$par

      }


      pop[now_year, 'effort'] <-
        distribute_fleet(
          pop = pop %>% filter(year == y),
          prior_profits = pop$profits[pop$year  == (y - 1)],
          year = y,
          burn_year = burn_year,
          effort = effort[y],
          fleet = fleet,
          num_patches = num_patches,
          mpa = mpa
        )

      if (fleet$tech_rate > 0 & y > burn_year){
        q[y] <- q[y - 1] + pmax(0,rnorm(1,fleet$tech_rate * q[y - 1], fleet$tech_rate*fleet$q))
      }

      pop[now_year, 'f'] <-
        pop[now_year, 'effort'] * q[y]

      # grow and die -----

      pop[pop$year == (y + 1), 'numbers'] <-
        pop[now_year,] %>%
        group_by(patch) %>%
        mutate(numbers = grow_and_die(
          numbers = numbers,
          f = f,
          mpa = mpa,
          fish = fish,
          fleet = fleet
        )$survivors) %>%
        ungroup() %>%
        {
          .$numbers
        }


      pop[now_year, 'numbers_caught'] <-
        pop[now_year,] %>%
        group_by(patch) %>%
        mutate(numbers_caught = grow_and_die(
          numbers = numbers,
          f = f,
          mpa = mpa,
          fish = fish,
          fleet = fleet
        )$caught) %>%

        ungroup() %>%
        {
          .$numbers_caught
        }
      # if (y > burn_year){ browser()}
      pop <- pop %>%
        # group_by(patch,year) %>%
        mutate(patch_age_costs = ((cost) * (effort) ^ fleet$beta) / fish$max_age) %>% #divide costs up among each age class
        # ungroup() %>%
        mutate(
          ssb = numbers * ssb_at_age,
          biomass = numbers * weight_at_age,
          biomass_caught = numbers_caught * weight_at_age,
          profits = biomass_caught * fish$price - patch_age_costs
        )

      # Adjust fleet

      if (y > burn_year & fleet$fleet_model == 'open-access') {
        effort[y + 1] <-
          max(0, effort[y] + fleet$theta * sum(pop$profits[now_year]))
      } else if (y > burn_year &
                 fleet$fleet_model == 'constant-effort') {
        effort[y + 1] <- effort[y]

      }

      # spawn ----

      pop$numbers[pop$year == (y + 1) &
                    pop$age == fish$min_age] <-
        calculate_recruits(
          pop = pop[pop$year == (y + 1), ],
          fish = fish,
          num_patches = num_patches,
          phase = model_phase,
          move_matrix = larval_move_matrix,
          rec_devs = rec_devs[y + 1]
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

        effort[y + 1] <- fleet$initial_effort

      }


    }

    pop <- pop %>%
      filter(year > burn_year, year < max(year)) %>%
      mutate(eventual_mpa = patch %in% mpa_locations)

    return(pop)

  }