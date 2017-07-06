#' \code{run_catch_curve} runs a weighted catch curve on spasm data
#'
#' @param length_comps a data frame with rows length_bin and numbers
#' @param fish a spasm fish object
#'
#' @return total mortality estimated by catch curve, z
#' @export
#'
#' @examples run_catch_curve(length_comps, fish)
run_catch_curve <- function(length_comps, fish) {

  age_comps <-  length_to_age(length_samples = length_comps,
                                 cv = fish$cv_len,
                                 k = fish$vbk,
                                 linf = fish$linf,
                                 t0 = fish$t0,
                                 max_age = fish$max_age)

  cc_dat <- age_comps %>%
    ungroup() %>%
    mutate(log_numbers = log(pmax(1e-3,numbers)))

  peak_age <- cc_dat$age[cc_dat$numbers == max(cc_dat$numbers)]

  cc_dat <- cc_dat %>%
    filter(age > peak_age)

  cc <- lm(log_numbers ~ age, data = cc_dat)

  ln_hat <- predict(cc) %>% as.numeric()

  cc_weights <- (ln_hat - min(ln_hat)) / sum(ln_hat - min(ln_hat))

  cc <- lm(log_numbers ~ age, data = cc_dat, weights = cc_weights)

  z <- cc$coefficients['age'] %>% as.numeric()

  return(z)

}