#' Title
#'
#' @param effort
#' @param fish
#' @param fleet
#' @param sim_years
#' @param burn_year
#' @param mpa_size
#' @param mpa_year
#' @param num_patches
#' @param use
#'
#' @return an estiamte of MSY
#' @export
#'
estimate_msy <-
  function(effort,
           fish,
           fleet,
           sim_years = 25,
           burn_year = 25,
           mpa_size = 0,
           mpa_year = 100,
           num_patches = 1,
           use = "fit") {

     fleet <-
      spasm::update_fleet(
        fleet = purrr::list_modify(
          fleet,
          fleet_model = "constant-effort",
          initial_effort = effort,
          sigma_effort = 0
        ),
        fish = fish
      )

    set.seed(24)

    sim <- spasm::sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(mpa_size = 0),
      num_patches = num_patches,
      sim_years = sim_years,
      burn_year = burn_year,
      time_step = fish$time_step,
      tune_costs = F,
      est_msy = F
    )

    yields <- sim %>%
      filter(year == max(year))

    if (use == "fit"){
      out <- -(sum(yields$biomass_caught))
    } else{

      out <- sum(yields$biomass)

    }

    return(out)

  }