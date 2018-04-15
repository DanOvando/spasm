#' Title
#'
#' @param cost
#' @param fish
#' @param fleet
#' @param msy
#' @param e_msy
#' @param b_msy
#' @param b_v_bmsy_target
#' @param p_response
#' @param sim_years
#' @param burn_year
#' @param num_patches
#' @param use
#'
#' @return an estimate of costs
#' @export
#'
estimate_costs <-
  function(cost,
           fish,
           fleet,
           msy,
           e_msy,
           b_msy,
           b_v_bmsy_oa,
           p_response,
           sim_years = 100,
           burn_year = 25,
           num_patches = 1,
           use = "fit") {


    fish$msy <- msy

    fleet$cost <- cost

    fleet$theta <- p_response

    fleet$e_msy <- e_msy

    fleet$p_msy <-
      fish$price * fish$msy - fleet$cost * fleet$e_msy ^ fleet$beta

    fleet$fleet_model = "open-access"
    fleet$initial_effort = 0.1 * (fish$m * 0.8)  / fleet$q
    fleet$sigma_effort = 0

    # fleet <-
    #   spasm::update_fleet(
    #     fleet = purrr::list_modify(
    #       fleet,
    #       fleet_model = "open-access",
    #       initial_effort = 10,
    #       sigma_effort = 0
    #     ),
    #     fish = fish
    #   )
    sim <- spasm::sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(mpa_size = 0),
      num_patches = num_patches,
      sim_years = sim_years,
      burn_year = burn_year,
      time_step = fish$time_step,
      est_msy = F,
      tune_costs = F
    )

    final_b <- sim %>%
      filter(year == max(year)) %>% {
        sum(.$biomass)
      }

    final_e <- sim %>%
      filter(year == max(year)) %>% {
        unique(.$effort)
      }

    print(final_b / b_msy)

    if (use == "fit") {
      out <- ((final_b / b_msy) - b_v_bmsy_oa) ^ 2
    } else{
      out <- final_ef
    }

  }