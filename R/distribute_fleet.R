#' \code{distribute_fleet} moves the fishing fleet around
#'
#' @param pop
#' @param effort
#' @param fleet
#'
#' @return a redistributed vector of effort in space
#' @export
#'
#' @examples distribute_fleet(pop, effort, fleet)
distribute_fleet <-
  function(pop,
           effort,
           prior_profits,
           year,
           burn_year,
           fleet,
           num_patches,
           mpa) {
    if (fleet$effort_allocation == 'simple') {
      pop$effort[pop$mpa == F] <-
        effort / length(unique(pop$patch[pop$mpa == F]))

      # pop$effort <- effort / num_patches

    }
    if (fleet$effort_allocation == 'gravity') {
      pop$effort[pop$mpa == F] <-
        effort * ((pop$numbers[pop$mpa == F] * pop$ssb_at_age[pop$mpa == F]) / sum((pop$numbers[pop$mpa == F] * pop$ssb_at_age[pop$mpa == F])))

    }

    if (fleet$effort_allocation == 'profit-gravity') {
      if (all(is.na(prior_profits[pop$mpa == F])) | year <= burn_year | all(prior_profits[pop$mpa == F] == 0))
      {
        pop$effort[pop$mpa == F] <-
          effort * ((pop$numbers[pop$mpa == F] * pop$ssb_at_age[pop$mpa == F]) / sum((pop$numbers[pop$mpa == F] * pop$ssb_at_age[pop$mpa == F])))
      } else {
        pop$effort[pop$mpa == F] <-
          effort * ((prior_profits[pop$mpa == F])^2 / sum((prior_profits[pop$mpa == F])^2, na.rm = T))

      }
    }

    return(pop$effort)


  }