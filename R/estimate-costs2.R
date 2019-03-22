#' estimate costs that produce a given B/Bmsy open access
#'
#' @param fish
#' @param fleet
#' @param msy
#' @param e_msy
#' @param b_msy
#' @param p_response
#' @param sim_years
#' @param burn_year
#' @param num_patches
#' @param use
#' @param max_cr_ratio the maximum cost to revenue ratio
#' @param b_v_bmsy_oa
#' @param lags
#'
#' @return an estimate of costs
#' @export
#'
estimate_costs2 <-
  function(pars,
           fish,
           fleet,
           b_ref_oa,
           sim_years = 100,
           burn_years = 25,
           num_patches = 1,
           use = "fit",
           lags = 0,
           sprinkler = FALSE,
           mpa_habfactor = 1) {

    fleet$max_cr_ratio <- pars[1]

    fleet$max_perc_change_f <- pars[2]


    # fleet$p_msy <-
    #   fish$price * fish$msy - fleet$cost * fleet$e_msy ^ fleet$beta

    # if (fleet$p_msy < 0){
    #   # print("d'oh")
    #
    #   if (use == "fit") {
    #     out <- ((fleet$p_msy) - b_v_bmsy_oa) ^ 2
    #   } else{
    #     out <- NA
    #   }
    #
    # } else {
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
      burn_years = burn_years,
      time_step = fish$time_step,
      est_msy = FALSE,
      tune_costs = FALSE,
      sprinkler = sprinkler,
      mpa_habfactor = mpa_habfactor
    )
    b0 <- sum(sim$biomass[sim$year == (burn_years + 1)])

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

    penalty <- sum(sim$effort[sim$year > burn_years] < 1)

    # write(final_b / b_msy, "wtf.txt", append = T)
    if (use == "fit") {
      out <- (log(final_b / b0) - log(b_ref_oa)) ^ 2 + penalty
    } else{
      out <- final_e
    }

    # } # close negative pmsy option

    return(out)

  }