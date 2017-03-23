#' create_fleet
#'
#' @param eq_f
#' @param length_50_sel
#' @param length_95_sel
#' @param mpa_reaction
#' @param fish
#'
#' @return
#' @export
#'
#' @examples create_fleet(eq_f = 2,length_50_sel = 25, length_95_sel = 27, fish = bluefish)
#'
create_fleet <- function(eq_f = NA,
                         length_50_sel = NA,
                         length_95_sel = NA,
                         mpa_reaction = 'concentrate',
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
    map_reaction = mpa_reaction
  )
}