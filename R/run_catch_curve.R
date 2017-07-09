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
  age_comps <-  length_to_age(
    length_samples = length_comps,
    cv = fish$cv_len,
    k = fish$vbk,
    linf = fish$linf,
    t0 = fish$t0,
    max_age = fish$max_age,
    min_age = fish$min_age,
    time_step = fish$time_step
  )

  cc_dat <- age_comps %>%
    ungroup() %>%
    mutate(log_numbers = ifelse(numbers > 0, log(numbers), NA))

  peak_age <- cc_dat$age[cc_dat$numbers == max(cc_dat$numbers)]

  # first_zero <- cc_dat$age[cc_dat$age > peak_age & cc_dat$numbers <= 1][1]
  #
  # if(is.na(first_zero)){first_zero <-  max(cc_dat$age) + 1}

  cc_dat <- cc_dat %>%
    filter(age >= peak_age)

  cc <- lm(log_numbers ~ age, data = cc_dat)

  pos_ages <- cc_dat$numbers > 0

  cc_weights <- rep(0, length(pos_ages))

  ln_hat <- predict(cc) %>% as.numeric()

  ln_hat <- (ln_hat - min(ln_hat)) / sum(ln_hat - min(ln_hat))

  cc_weights[pos_ages] <- ln_hat

  cc <- lm(log_numbers ~ age, data = cc_dat, weights = cc_weights)

  z <- (cc$coefficients['age'] %>% as.numeric()) * fish$time_step

  return(z)

}