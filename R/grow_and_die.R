#' \code{grow_and_die} runs growth mortality and fishing
#'
#' @param numbers
#' @param fish
#' @param fleet
#'
#' @return numbers and catch at age for all ages above recruits
#' @export
#'
#' @examples grow_and_die(numbers, fish, fleet)
grow_and_die <- function(numbers, f, mpa, fish, fleet) {
  survivors <- vector(mode = 'numeric', length = length(numbers))

  # survival  <- exp(-(fish$m + (f * (!mpa) * fleet$sel_at_age)))

  survival  <- exp(-fish$time_step*(fish$m + (f * fleet$sel_at_age)))

  death <-  1 - survival

  max_index <- length(survivors)

  survivors[2:max_index] <-
    numbers[1:(max_index - 1)] * survival[1:(max_index - 1)]

  survivors[max_index] <-
    survivors[max_index] + numbers[max_index] * survival[max_index]

  caught <-
    (fish$time_step * f * fleet$sel_at_age) / (fish$time_step * (fish$m + (f * fleet$sel_at_age))) *  (numbers * death)
  # return(survivors)

  return(list(survivors = survivors, caught = caught))


}
