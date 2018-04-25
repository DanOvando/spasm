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
           use = "fit",
           lags = 0) {


    fish$msy <- msy

    fleet$cost <- cost

    fleet$theta <- p_response

    fleet$e_msy <- e_msy

    fleet$p_msy <-
      fish$price * fish$msy - fleet$cost * fleet$e_msy ^ fleet$beta

    if (fleet$p_msy < 0){
      # print("d'oh")

      if (use == "fit") {
        out <- ((fleet$p_msy) - b_v_bmsy_oa) ^ 2
      } else{
        out <- NA
      }

    } else {
    fleet$fleet_model = "open-access"
    # fleet$initial_effort = 0.1 * (fish$m * 0.8)  / fleet$q
    fleet$sigma_effort = 0

    set.seed(42)
    sim <- spasm::sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(mpa_size = 0),
      num_patches = num_patches,
      sim_years = sim_years,
      burn_year = burn_year,
      time_step = fish$time_step,
      est_msy = F,
      tune_costs = F,
      b_v_bmsy_oa = b_v_bmsy_oa
    )

    final_b <- sim %>%
      filter(year >= (max(year) - lags)) %>%
      group_by(year) %>%
      summarise(biomass = sum(biomass)) %>%
      ungroup() %>%
      {
        mean(.$biomass)
      }

    final_e <- sim %>%
      filter(year >= (max(year) - lags)) %>%
      group_by(year) %>%
      summarise(effort = mean(effort)) %>%
      ungroup() %>%
      {
        mean(.$effort)
      }

    # write(final_b / b_msy, "wtf.txt", append = T)

    if (use == "fit") {
      out <- ((final_b / b_msy) - b_v_bmsy_oa) ^ 2
    } else{
      out <- final_ef
    }

    } # close negative pmsy option

    return(out)

  }