#' create_fleet
#'
#' @param eq_f
#' @param length_50_sel
#' @param length_95_sel
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
#'
#' @return a fleet object
#' @export
#'
#' @examples create_fleet(eq_f = 2,length_50_sel = 25, length_95_sel = 27, fish = bluefish)
#'
create_fleet <- function(eq_f = NA,
                         length_50_sel = 1,
                         length_95_sel = 2,
                         mpa_reaction = 'concentrate',
                         price = 1,
                         cost = .1,
                         beta = 1.3,
                         theta = 1e-1,
                         q = 1e-3,
                         fleet_model = 'constant-effort',
                         effort_allocation = 'gravity',
                         initial_effort = 100,
                         target_catch = 0,
                         fish) {


  age_50_sel <- (log(1 - length_50_sel / fish$linf) / -fish$vbk) + fish$t0

  age_95_sel <- (log(1 - length_95_sel / fish$linf) / -fish$vbk) + fish$t0

  sel_at_age <-
  ((1 / (1 + exp(-log(
    19
  ) * (((1:fish$max_age) - age_50_sel) / (age_95_sel - age_50_sel)
  )))))

  fleet <- list(
    eq_f = eq_f,
    length_50_sel = length_50_sel,
    length_95_sel = length_95_sel,
    sel_at_age = sel_at_age,
    mpa_reaction = mpa_reaction,
    price = price,
    cost = cost,
    beta = beta,
    theta = theta,
    q = q,
    fleet_model = fleet_model,
    effort_allocation = effort_allocation,
    initial_effort = initial_effort,
    target_catch = target_catch
  )
}