#' Generate a distribution of values over time
#'
#' Generates a time series of values with lognormal errors and
#' some temporal autocorrelation
#'
#' @param thing the thing to be jittered
#' @param sigma the degree of variation
#' @param ac the degree of autocorrelation 0-1
#' @param time the length of time series to simulate
#'
#' @return a vector of length time
#' @export
#'
#' @examples
#' generate_timeseries(thing = 2, sigma = .1, ac = 0.5, time = 10)
generate_timeseries <- function(thing, sigma, ac, time){

  if (length(thing) == 1 & sigma > 0){

    thing <- thing * exp(rnorm(time, 0, sigma) - sigma^2/2)

    for (t in 2:length(thing)) {

      thing[t] <-
        thing[t - 1] * ac + sqrt(1 - ac ^ 2) * thing[t]
    }

  }

  return(thing)

}