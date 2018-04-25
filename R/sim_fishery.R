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
           tune_costs = F,
           b_v_bmsy_oa = 0.5,
           time_step,
           max_window = 10) {

    if (est_msy == T){

      tol <- 100

      lower <- 0

      upper <- 400

      golden <- (sqrt(5) -1)/2

      best <- 1000

      delta_best <- 100

      counter <-  0
#
#       efforts <- seq(0,10000, by = 50)
#
#       huh <- NA
#
#       for (i in seq_along(efforts)){
#
#         huh[i] <- estimate_msy(efforts[i], fish = fish, fleet = fleet)
#
#       }
#       # browser()

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

        # if (counter == 20){
        #   browser()
        # }

      } # close golden while

     msy_fit <- nlminb(mean(c(lower, upper)), estimate_msy, fish = fish, fleet = fleet, lower = 0)

      fleet$e_msy <- msy_fit$par

      fish$msy <- -msy_fit$objective

      fish$b_msy <- estimate_msy(fleet$e_msy, fish = fish, fleet = fleet, use = "other")

      msy <- fish$msy

    }

    if (tune_costs == T){
      tol <- .01

      lower <- 0

      upper <- 100

      golden <- (sqrt(5) -1)/2

      best <- 1000

      delta_best <- 100

      counter <- 0

      set.seed(24)

      while(delta_best > tol) {

        counter <- counter + 1

        constant <- (1 - golden) * (upper - lower)

        x1 <- lower + constant

        x2 <- upper - constant

        ss_1 <- estimate_costs(cost = x1,
                                fish = fish,
                                fleet = fleet,
                                msy = fish$msy,
                                e_msy = fleet$e_msy,
                                b_msy = fish$b_msy,
                                p_response = fleet$theta,
                                b_v_bmsy_oa = b_v_bmsy_oa,
                               sim_years = 100)

        ss_2 <- estimate_costs(cost = x2,
                               fish = fish,
                               fleet = fleet,
                               msy = fish$msy,
                               e_msy = fleet$e_msy,
                               b_msy = fish$b_msy,
                               p_response = fleet$theta,
                               b_v_bmsy_oa = b_v_bmsy_oa,
                               sim_years = 100)

        delta_best <-  (0 -  min(ss_1,ss_2))^2

        best <- min(ss_1,ss_2)

        if (ss_1 < ss_2){

          lower <- lower

          upper <- x2
        } else{

          lower <- x1

          upper <- upper

        }

        if (counter > 20){
          delta_best <- 0
        }

      } # close golden while

cost_fit <-         nlminb(
  mean(c(lower, upper)),
  estimate_costs,
  fish = fish,
  fleet = fleet,
  lower = 0,
  msy = fish$msy,
  e_msy = fleet$e_msy,
  b_msy = fish$b_msy,
  p_response = fleet$theta,
  b_v_bmsy_oa = b_v_bmsy_oa
)
      fleet$cost <- cost_fit$par


    } # close estimate costs

    fleet$p_msy <- fish$price * fish$msy - fleet$cost * fleet$e_msy ^ fleet$beta

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
      rec_devs <-
        rnorm(
          sim_years,
          mean = -(fish$sigma_r ^ 2) / 2,
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
      fish$r0 / num_patches * exp(-fish$m * seq(fish$min_age, fish$max_age, fish$time_step))

    n0_at_age[fish$max_age] <-
      n0_at_age[fish$max_age] / (1 - exp(-fish$m))

    b0_at_age <- n0_at_age * fish$weight_at_age

    ssb0_at_age <- n0_at_age * fish$ssb_at_age


    # generate time series of price, cost, and q if called for
    price <- generate_timeseries(fish$price, sigma = fish$price_cv * fish$price, ac = fish$price_ac, time = sim_years)

    q <- generate_timeseries(fleet$q, sigma = fleet$q_cv * fleet$q, ac = fleet$q_ac, time = sim_years)

    #
#     if (tune_costs == T){
#
#       fleet$cost <- fleet$oa_ratio * sum(ssb0_at_age) * mean(price) * mean(q)
#       cost <- generate_timeseries(fleet$cost, sigma = fleet$cost_cv * fleet$cost, ac = fleet$cost_ac, time = sim_years)
#
#     } else{

    cost <- generate_timeseries(fleet$cost, sigma = fleet$cost_cv * fleet$cost, ac = fleet$cost_ac, time = sim_years)

    # }

    if (length(q) == 1) {
      q <- rep(q, sim_years)
    }
    if (length(price) == 1) {
      price <- rep(price, sim_years)
    }

    cost_frame <- data_frame(year = 1:sim_years, cost = cost)

    pop <- pop %>%
      select(-cost) %>%
      left_join(cost_frame, by = "year")

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

      pop[now_year &
        pop$age > fish$min_age, ] <-
        move_fish(
          pop %>% filter(year == y, age > fish$min_age),
          fish = fish,
          num_patches = num_patches,
          move_matrix = adult_move_matrix
        )

      # change management

      if ((y - burn_year) == manager$year_mpa) {
        prop_mpas <- floor(num_patches * manager$mpa_size)

        if (random_mpas == T) {
          mpa_locations <- sample(1:num_patches, prop_mpas)
        } else {
          mpa_locations <- (1:num_patches)[0:prop_mpas]
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

          # new_theta <- (effort[y - 1] * fleet$theta_tuner) / (total_initial_profits + 1e-3)
          #
          # if (new_theta < 0) {
          #   stop("fishery is unprofitable at b0")
          # }

          # fleet <- purrr::list_modify(fleet, theta = new_theta)

          # fleet <- update_fleet(fleet = purrr::list_modify(fleet, theta = new_theta), fish = fish)
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
          e_msy = fleet$e_msy,
          p_msy = fleet$p_msy,
          mey_buffer = fleet$mey_buffer,
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
      # if (y > burn_year){ browser()}
      pop <- pop %>%
        # group_by(patch,year) %>%
        mutate(patch_age_costs = ((cost) * (effort) ^ fleet$beta) / fish$max_age) %>% # divide costs up among each age class
        # ungroup() %>%
        mutate(
          ssb = numbers * ssb_at_age,
          biomass = numbers * weight_at_age,
          biomass_caught = numbers_caught * weight_at_age,
          profits = biomass_caught * price[y] - patch_age_costs
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

    price_mat <- data_frame(year = 1:sim_years, price = price)

    pop <- pop %>%
      left_join(rec_mat, by = "year") %>%
      left_join(enviro_mat, by = "year") %>%
      left_join(price_mat, by = "year") %>%
      filter(year > burn_year, year < max(year)) %>%
      mutate(eventual_mpa = patch %in% mpa_locations)

    return(pop)
  }