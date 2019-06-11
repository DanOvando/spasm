#' determine_effort
#'
#' \code{determine_effort} determines the effort in the current time period based
#' on different functional forms
#'
#' @param last_effort effort in the last time period
#' @param fleet the fleet object
#' @param fish the fish object
#' @param y the current year
#' @param burn_years the last year of the burn period
#' @param pop the population object
#' @param mpa the mpa object
#' @param num_patches the number of patches in the system
#'
#' @return the effort in the current year
#' @export
#'
determine_effort <-
  function(last_effort,
           fleet,
           fish,
           y,
           burn_years,
           pop,
           mpa,
           num_patches,
           effort_devs,
           profit_lags = 4,
           previous_max = NA) {
    new_effort <- last_effort
    if (fleet$fleet_model == 'constant-catch') {
      effort_for_catch <- nlminb(
        1,
        catch_target,
        target_catch =
          fleet$target_catch,
        pop = pop %>% filter(year == y),
        num_patches = num_patches,
        mpa = mpa,
        fleet = fleet,
        lower = 0,
        use = 'opt',
        fish = fish,
        prior_profits = pop$profits[pop$year == (y - 1)],
        year = y,
        burn_years = burn_years
      )

      new_effort <- effort_for_catch$par


    }

    if (fleet$fleet_model == 'supplied-catch') {
      target_catch <- fleet$catches[y - burn_years]

      effort_for_catch <- nlminb(
        1,
        catch_target,
        target_catch =
          target_catch,
        pop = pop %>% filter(year == y),
        num_patches = num_patches,
        mpa = mpa,
        fleet = fleet,
        lower = 0,
        use = 'opt',
        fish = fish,
        prior_profits = pop$profits[pop$year == (y - 1)],
        year = y,
        burn_years = burn_years
      )
      new_effort <- effort_for_catch$par

    }


    if (fleet$fleet_model == 'open-access') {
      profits <- pop %>%
        filter(year >= (y - (1 + profit_lags)), year < y) %>%
        group_by(year) %>%
        summarise(profits = sum(profits),
                  effort = sum(effort[age == 0])) %>%
        mutate(ppue = profits / (effort + 1e-6))

      # new_effort <-
      #   pmin(
      #     last_effort * (1 + fleet$max_perc_change_f),
      #     last_effort + (fleet$theta * weighted.mean(profits$ppue, profits$year)) * exp(effort_devs[y + 1])
      #   )
      new_effort <-
          last_effort + (fleet$theta * weighted.mean(profits$ppue, profits$year)) * exp(effort_devs[y + 1])

      new_effort[is.na(new_effort)] <- 1e-3

      if (new_effort <= 1e-3) {
        new_effort = 1e-3 / (2 - new_effort / 1e-3)

      }

    }

    if (fleet$fleet_model == 'constant-effort') {
      new_effort <- last_effort

    }

    if (fleet$fleet_model != 'open-access') {
      new_effort <- new_effort * exp(effort_devs[y + 1])
    }
    if (fleet$fleet_model == "random_walk") {
      new_effort <- last_effort * exp(effort_devs[y + 1])

    }

    # new_effort <- pmax(last_effort * (1 - fleet$max_perc_change_f), pmin(new_effort, previous_max * (1 + fleet$max_perc_change_f)))

    return(new_effort)


  }