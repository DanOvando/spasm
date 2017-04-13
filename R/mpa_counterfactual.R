#' \code{mpa_counterfactual} compares a fishery with and without an MPA
#'
#' @param fish
#' @param fleet
#' @param year_mpa
#' @param mpa_size
#' @param sim_years
#' @param num_patches
#' @param burn_years
#'
#' @return an object comparing runs
#' @export
#'
mpa_counterfactual <- function(fish, fleet, year_mpa,mpa_size,sim_years, num_patches,
                     burn_years){


  no_mpa <- sim_fishery(fish = fish, fleet = fleet, manager = create_manager(year_mpa = year_mpa, mpa_size = 0), sim_years = sim_years, num_patches = num_patches,
                        burn_years = burn_years) %>%
    mutate(experiment = 'no-mpa')

  wi_mpa <- sim_fishery(fish = fish, fleet = fleet, manager = create_manager(year_mpa = year_mpa, mpa_size = mpa_size), sim_years = sim_years, num_patches = num_patches,
                        burn_years = burn_years) %>%
    mutate(experiment = 'with-mpa')

  outcomes <- no_mpa %>%
    bind_rows(wi_mpa) %>%
    group_by(year,experiment) %>%
    summarise(ssb = sum(ssb), percent_mpa = mean(mpa),
              catch = sum(biomass_caught),
              profits = sum(profits),
              effort = sum(effort))

  return(list(outcomes = outcomes, wi_mpa = wi_mpa))



}