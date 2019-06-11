#' \code{distribute_fleet} moves the fishing fleet around
#'
#' @param pop
#' @param effort
#' @param fleet
#'
#' @return a redistributed vector of effort in space
#' @export
#'
#' @examples
#' \dontrun{
#'   distribute_fleet(pop, effort, fleet)
#'   }
distribute_fleet <-
  function(pop,
           effort,
           prior_profits,
           year,
           burn_years,
           fleet,
           num_patches,
           mpa) {
    if (fleet$effort_allocation == 'simple') {
      pop$effort[pop$mpa == F] <-
        effort / length(unique(pop$patch[pop$mpa == F]))

      # pop$effort <- effort / num_patches

    }
    if (fleet$effort_allocation == 'gravity') {
      ssb_by_patch <-  pop %>%
        group_by(patch) %>%
        summarise(mpa = unique(mpa),
                  ssb = sum(numbers * ssb_at_age))

      effort_by_patch <-
        effort * (ssb_by_patch$ssb[ssb_by_patch$mpa == F] / sum(ssb_by_patch$ssb[ssb_by_patch$mpa == F])) %>%
        rep(each = length(unique(pop$age)))

      pop$effort[pop$mpa == F] <-  effort_by_patch

    }

    if (fleet$effort_allocation == 'profit-gravity') {
      if (all(is.na(prior_profits[pop$mpa == F])) |
          year <= burn_years | all(prior_profits[pop$mpa == F] == 0))
      {

         ssb_by_patch <-  pop %>%
          group_by(patch) %>%
           filter(mpa == F) %>%
          summarise(ssb = sum(numbers * ssb_at_age))

        effort_by_patch <- effort * (ssb_by_patch$ssb / sum(ssb_by_patch$ssb)) %>%
          rep(each = length(unique(pop$age)))

        pop$effort[pop$mpa == F] <-  effort_by_patch

              } else {

                profits_by_patch <-  pop %>%
                  select(patch,mpa) %>%
                  bind_cols(prior_profits = prior_profits) %>%
                  mutate(prior_profits = prior_profits - (min(prior_profits) - 1e-3)) %>%
                  group_by(patch) %>%
                  filter(mpa == F) %>%
                  summarise(profits = sum(prior_profits))

                effort_by_patch <- effort * (profits_by_patch$profits / sum(profits_by_patch$profits)) %>%
                  rep(each = length(unique(pop$age)))

                pop$effort[pop$mpa == F] <-  effort_by_patch

      }
    }

    return(pop$effort)


  }