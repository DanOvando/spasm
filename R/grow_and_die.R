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

  survival  <- exp(-(fish$m + (f * fleet$sel_at_age)))

  death <-  1 - survival

  survivors[2:fish$max_age] <-
    numbers[1:(fish$max_age - 1)] * survival[1:(fish$max_age - 1)]

  survivors[fish$max_age] <-
    survivors[fish$max_age] + numbers[fish$max_age] * survival[fish$max_age]

  caught <-
    (f * fleet$sel_at_age) / (fish$m + (f * fleet$sel_at_age)) *  (numbers * death)
  # return(survivors)

  return(list(survivors = survivors, caught = caught))


}
