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
           random_mpas = F,
           enviro = NA,
           enviro_strength = 1,
           rec_driver = "stochastic",
           est_msy = F,
           time_step,
           max_window = 10) {

    msy <- NA

    p_msy <- NA

    e_msy <- NA

    max_r_msy <-  NA

    if (est_msy == T){


      e_msy_guess <- fish$m/fleet$q

      tol <- 100

      lower <- 0

      upper <- e_msy_guess * 3

      golden <- (sqrt(5) -1)/2

      best <- 1000

      delta_best <- 100

      counter <-  0

      while(delta_best > tol | counter < 20) {

        counter <- counter + 1

        constant <- (1 - golden) * (upper - lower)

        x1 <- lower + constant

        x2 <- upper - constant

        yield_1 <- estimate_msy(x1, fish = fish, fleet = fleet)

        yield_2 <- estimate_msy(x2, fish = fish, fleet = fleet)

        delta_best <-  (best -  min(yield_1,yield_2))^2

        best <- min(yield_1,yield_2)

        if (yield_1 < yield_2){

          lower <- lower

          upper <- x2
        } else{

          lower <- x1

          upper <- upper

        }


      } # close golden while


      msy_foo <- function(seed, lower, upper, fish, fleet) {

        msy_fit <-
          nlminb(
            mean(c(lower, upper)),
            estimate_msy,
            fish = fish,
            fleet = fleet,
            lower = 0,
            seed = seed
          )

        out <- list()

        out$e_msy <- msy_fit$par

        out$msy <- -msy_fit$objective

        ref_points <-
          estimate_msy(
            out$e_msy,
            fish = fish,
            fleet = fleet,
            use = "other",
            seed = seed
          )

        out$b_msy <- ref_points$b_msy

        out$r_msy <- ref_points$r_msy

        return(out)

      } # close msy foo

      msy_ests <- map(sample(1:10000, 10, replace = F), msy_foo, fish = fish, fleet = fleet, lower = lower, upper = upper)

      max_r_msy <- max(map_dbl(msy_ests,"r_msy"))

      e_msy <- mean(map_dbl(msy_ests,"e_msy"))

      fleet$e_msy <- e_msy

      msy <- mean(map_dbl(msy_ests,"msy"))

      fish$msy <- msy

    }

    sim_years <- burn_year + sim_years

    if (fleet$fleet_model == "supplied-catch") {
      sim_years <- sim_years + 1
    }


    pop <-
      expand.grid(
        year = 1:sim_years,
        patch = 1:num_patches,
        age = seq(fish$min_age, fish$max_age, fish$time_step)
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
        cost = NA
      ) %>%
      as_data_frame() %>%
      arrange(year, patch, age)


    effort <- vector(mode = "double", length = sim_years)

    f <- vector(mode = "double", length = sim_years)


    if (rec_driver == "stochastic") {
      # rec_devs <-
      #   rnorm(
      #     sim_years,
      #     mean = -(fish$sigma_r ^ 2) / 2,
      #     sd = fish$sigma_r
      #   )

      rec_devs <-
        rnorm(
          sim_years,
          mean = 0,
          sd = fish$sigma_r
        )

      ## autocorrelated recruitment deviations
      for (t in 2:length(rec_devs)) {
        rec_devs[t] <-
          rec_devs[t - 1] * fish$rec_ac + sqrt(1 - fish$rec_ac ^ 2) * rec_devs[t]
      }
    } else if (rec_driver == "environment") {
      if (length(enviro) != sim_years) {
        stop("environment must be same length as sim_years")
      }

      rec_devs <-
        rnorm(
          sim_years,
          mean = enviro_strength * enviro,
          sd = fish$sigma_r
        )
    }

    effort_devs <-
      rnorm(
        sim_years,
        mean = 0,
        sd = fleet$sigma_effort
      )

    for (t in 2:length(effort_devs)) {
      effort_devs[t] <-
        effort_devs[t - 1] * fleet$effort_ac + sqrt(1 - fleet$effort_ac ^ 2) * effort_devs[t]
    }

    mpa_locations <- -1

    n0_at_age <-
      (fish$r0 / num_patches) * exp(-fish$m * seq(fish$min_age, fish$max_age, fish$time_step))

    n0_at_age[fish$max_age + 1] <-
      n0_at_age[fish$max_age + 1] / (1 - exp(-fish$m))

    b0_at_age <- n0_at_age * fish$weight_at_age

    ssb0_at_age <- n0_at_age * fish$ssb_at_age

    # generate time series of price, cost, and q if called for

    price_series <- generate_timeseries(fish$price, cv = fish$price_cv, ac = fish$price_ac, percent_slope = fish$price_slope, time = sim_years)

    q <- generate_timeseries(fleet$q, cv = fleet$q_cv, ac = fleet$q_ac, percent_slope = fleet$q_slope, time = sim_years)

    if (length(q) == 1) {
      q <- rep(q, sim_years)
    }
    if (length(price_series) == 1) {
      price_series <- rep(price_series, sim_years)
    }

  # tune costs based on some heavy fishing at b0

   hyp_f <- fish$m #hypothetical f

   hyp_effort <- hyp_f / mean(q[(burn_year + 1):sim_years])

   hyp_f_at_age <- hyp_f * fleet$sel_at_age

   hyp_b0_catch <- sum((hyp_f_at_age / (hyp_f_at_age + fish$m))  * b0_at_age * (1 - exp(-(hyp_f_at_age + fish$m))))

   b0_revenue <- mean(price_series[(burn_year + 1):sim_years]) * hyp_b0_catch

   hyp_profits_guess <- b0_revenue * (1 - fleet$max_cp_ratio)

   cost_guess <-
     (b0_revenue - hyp_profits_guess) / hyp_effort ^ fleet$beta


   fleet$theta <- (fleet$max_perc_change_f * hyp_effort) / (hyp_profits_guess / hyp_effort)

    cost_series <- generate_timeseries(cost_guess, cv = fleet$cost_cv, ac = fleet$cost_ac, percent_slope = fleet$cost_slope, time = sim_years)

    cost_series <- (cost_series / max(cost_series)) * cost_guess

    fleet$cost <- cost_guess

    price_frame <- data_frame(year = 1:sim_years, price = price_series)

    cost_frame <- data_frame(year = 1:sim_years, cost = cost_series)

    pop <- pop %>%
      select(-cost) %>%
      left_join(cost_frame, by = "year") %>%
      left_join(price_frame, by = "year")

    pop$numbers[pop$year == 1] <- rep(n0_at_age, num_patches)

    if (fleet$cost_function == "distance from port") {
      cost_frame <- expand.grid(year = 1:sim_years, patch = 1:num_patches) %>%
        as_data_frame() %>%
        left_join(cost_frame, by = "year") %>%
        mutate(cost = cost * (1 + fleet$cost_slope * (patch - 1)))

      pop <- pop %>%
        select(-cost) %>%
        left_join(cost_frame, by = c("patch", "year"))
    }


    pop <- pop %>%
      left_join(
        data_frame(
          age = seq(fish$min_age, fish$max_age, fish$time_step),
          ssb_at_age = fish$ssb_at_age,
          weight_at_age = fish$weight_at_age
        ),
        by = "age"
      )

    y <- 1

    model_phase <- "burn"

    adult_move_grid <-
      expand.grid(
        source = 1:num_patches,
        sink = 1:num_patches
      ) %>%
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
      expand.grid(
        source = 1:num_patches,
        sink = 1:num_patches
      ) %>%
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

      if (fish$density_movement_modifier > 0) {
      adult_density_modifier <- calc_density_gradient(pop, y, num_patches = num_patches,
                                                      density_modifier = fish$density_movement_modifier)

      if (any(adult_density_modifier < 0)){browser()}
      } else {
        adult_density_modifier <- 1
      }

      if (num_patches > 1) {
      pop[now_year &
        pop$age > fish$min_age, ] <-
        move_fish(
          pop %>% filter(year == y, age > fish$min_age),
          fish = fish,
          num_patches = num_patches,
          move_matrix = (adult_move_matrix * adult_density_modifier) / rowSums(adult_move_matrix * adult_density_modifier)
        )
      }

      # change management

      if ((y - burn_year) == manager$year_mpa) {
        prop_mpas <- floor(num_patches * manager$mpa_size)

        if (random_mpas == T) {
          mpa_locations <- sample(1:num_patches, prop_mpas)
        } else {
          mpa_locations <- (1:num_patches)[0:prop_mpas] #weird zero is in case prop_mpas is zero
        }

        pop$mpa[pop$patch %in% mpa_locations & pop$year >= y] <- T
      }
      # Adjust fleet
      if (y > (burn_year)) {
        if (y == (burn_year + 2) & fleet$fleet_model == "open-access") {

          profits <- pop %>%
            filter(year >= (y - (1 + fleet$profit_lags)), year < y) %>%
            group_by(year) %>%
            summarise(profits = sum(profits))

          total_initial_profits <- mean(profits$profits)

        }

        previous_max <- ifelse(y > (burn_year + 1),max(effort[max(1,(y - 1 - max_window)):(y - 1)]),fleet$initial_effort)

        effort[y] <- determine_effort(
          last_effort = ifelse(y > (burn_year + 1), effort[y - 1], fleet$initial_effort),
          fleet = fleet,
          fish = fish,
          y = y,
          burn_year = burn_year,
          pop = pop,
          mpa = mpa,
          num_patches = num_patches,
          effort_devs = effort_devs,
          profit_lags = fleet$profit_lags,
          previous_max = previous_max
        )

      }

      pop[now_year, "effort"] <-
        distribute_fleet(
          pop = pop %>% filter(year == y),
          prior_profits = pop$profits[pop$year == (y - 1)],
          year = y,
          burn_year = burn_year,
          effort = effort[y],
          fleet = fleet,
          num_patches = num_patches,
          mpa = mpa
        )

      if (fleet$tech_rate > 0 & y > burn_year) {
        q[y] <-
          q[y - 1] + pmax(
            0,
            rnorm(1, fleet$tech_rate * q[y - 1], fleet$tech_rate * fleet$q)
          )
      }

      pop[now_year, "f"] <-
        pop[now_year, "effort"] * q[y]

      # grow and die -----

      pop[pop$year == (y + 1), "numbers"] <-
        pop[now_year, ] %>%
        group_by(patch) %>%
        mutate(numbers = grow_and_die(
          numbers = numbers,
          f = f,
          mpa = mpa,
          fish = fish,
          fleet = fleet,
          y = y
        )$survivors) %>%
        ungroup() %>%
        {
          .$numbers
        }


      pop[now_year, "numbers_caught"] <-
        pop[now_year, ] %>%
        group_by(patch) %>%
        mutate(
          numbers_caught = grow_and_die(
            numbers = numbers,
            f = f,
            mpa = mpa,
            fish = fish,
            fleet = fleet,
            y = y
          )$caught
        ) %>%
        ungroup() %>%
        {
          .$numbers_caught
        }

      pop <- pop %>%
        mutate(patch_age_costs = ((cost) * (effort) ^ fleet$beta) / fish$max_age) %>% # divide costs up among each age class
        mutate(
          ssb = numbers * ssb_at_age,
          biomass = numbers * weight_at_age,
          biomass_caught = numbers_caught * weight_at_age,
          profits = biomass_caught * price - patch_age_costs
        )

      # spawn ----

      # if (is.na(spawning_season) | ((((year) - floor(year))/spawning_season) == 1))

      pop$numbers[pop$year == (y + 1) &
        pop$age == fish$min_age] <-
        calculate_recruits(
          pop = pop[pop$year == y, ],
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
          ungroup() %>%
          {
            (.$ssb)
          }

        model_phase <- "recruit"

        effort[y + 1] <- fleet$initial_effort
      }
    }
    rec_mat <- data_frame(year = 1:sim_years, rec_dev = rec_devs)

    enviro_mat <- data_frame(year = 1:sim_years, enviro = enviro)

    pop <- pop %>%
      left_join(rec_mat, by = "year") %>%
      left_join(enviro_mat, by = "year") %>%
      filter(year > burn_year, year < max(year)) %>%
      mutate(eventual_mpa = patch %in% mpa_locations,
             msy = msy,
             p_msy = p_msy,
             e_msy = e_msy)

    return(pop)
  }