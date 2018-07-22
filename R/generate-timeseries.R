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
generate_timeseries <- function(thing, cv, ac, time){

  # add random walk with drift https://rpubs.com/ericnovik/ar1stan

  sigma <- sqrt(log(cv^2 + 1))

  if (length(thing) == 1 & sigma > 0){

    thing_devs <-
      rnorm(
        time,
        mean = -(sigma ^ 2) / 2,
        sd = sigma
      )

    for (t in 2:length(thing_devs)) {

      thing_devs[t] <-
        thing_devs[t - 1] * ac + sqrt(1 - ac ^ 2) * thing_devs[t]
    }

    thing <- thing * exp(thing_devs)

  }

  return(thing)

}