#' create_fleet
#'
#' @param eq_f
#' @param length_50_sel
#' @param delta cm above length50 at 95 selectivity
#' @param mpa_reaction
#' @param fish
#' @param price
#' @param cost
#' @param beta
#' @param theta
#' @param q
#' @param fleet_model
#' @param effort_allocation
#' @param initial_effort
#' @param cost_function
#' @param cost_slope
#' @param tech_rate
#' @param target_catch
#' @param catches
#' @param sigma_effort
#'
#' @return a fleet object
#' @export
#'
#' @examples create_fleet(eq_f = 2,length_50_sel = 25, length_95_sel = 27, fish = bluefish)
create_fleet <- function(eq_f = NA,
                         length_50_sel = 1,
                         delta = 2,
                         fish,
                         mpa_reaction = 'concentrate',
                         price = 1,
                         cost = .1,
                         beta = 1.3,
                         theta = 1e-1,
                         q = 1e-3,
                         fleet_model = 'constant-effort',
                         effort_allocation = 'gravity',
                         cost_function = 'constant',
                         cost_slope = 0,
                         tech_rate = 0,
                         initial_effort = 100,
                         target_catch = 0,
                         catches = NA,
                         sigma_effort = 0) {

  p_selected <- function(mu, sigma, l50, delta){

    length_dist <- pmax(0,rnorm(1000, mu, sigma))

    sel_dist <-   ((1 / (1 + exp(-log(
      19
    ) * ((length_dist - l50) / (delta)
    )))))

    mean_p_selected <- mean(sel_dist)

  }

  sel_at_age <- data_frame(age = 0:fish$max_age) %>%
    mutate(mean_length_at_age = fish$length_at_age) %>%
    mutate(sd_at_age = mean_length_at_age * fish$cv_len) %>%
    mutate(mean_sel_at_age = map2_dbl(mean_length_at_age, sd_at_age, p_selected, l50 = length_50_sel, delta = delta))


  length_95_sel <- (length_50_sel + delta)


  fleet <- list(
    eq_f = eq_f,
    length_50_sel = length_50_sel,
    delta = delta,
    sel_at_age = sel_at_age$mean_sel_at_age,
    mpa_reaction = mpa_reaction,
    price = price,
    cost = cost,
    beta = beta,
    theta = theta,
    q = q,
    fleet_model = fleet_model,
    effort_allocation = effort_allocation,
    initial_effort = initial_effort,
    target_catch = target_catch,
    catches = catches,
    cost_function = cost_function,
    cost_slope = cost_slope,
    tech_rate = tech_rate,
    sigma_effort = sigma_effort
  )
}
