#' determine_effort
#'
#' \code{determine_effort} determines the effort in the current time period based
#' on different functional forms
#'
#' @param last_effort effort in the last time period
#' @param fleet the fleet object
#' @param fish the fish object
#' @param y the current year
#' @param burn_year the last year of the burn period
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
           burn_year,
           pop,
           mpa,
           num_patches,
           effort_devs,
           profit_lags = 4,
           e_msy,
           p_msy,
           mey_buffer = 2,
           previous_max = NA,
           max_expansion = 1.1) {
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
        burn_year = burn_year
      )

      new_effort <- effort_for_catch$par


    }

    if (fleet$fleet_model == 'supplied-catch') {
      target_catch <- fleet$catches[y - burn_year]

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
        burn_year = burn_year
      )
      new_effort <- effort_for_catch$par

    }


    if (fleet$fleet_model == 'open-access') {
      profits <- pop %>%
        filter(year >= (y - (1 + profit_lags)), year < y) %>%
        group_by(year) %>%
        summarise(profits = sum(profits))

      if (is.na(e_msy) | is.na(p_msy)) {
        stop("need to estiamte msy and tune costs to run open-access")
      }

      new_effort <-
        last_effort + e_msy * (fleet$theta * mean(profits$profits / (p_msy * mey_buffer))) * exp(effort_devs[y + 1])

      if (new_effort <= 0) {
        new_effort = -1 / (new_effort - 1)

      }

    }

    if (fleet$fleet_model == 'constant-effort') {
      new_effort <- last_effort

    }

    if (fleet$fleet_model != 'open-access'){
    new_effort <- new_effort * exp(effort_devs[y + 1])
    }
    new_effort <- pmin(new_effort, previous_max * max_expansion)

    return(new_effort)


  }