grow_and_die <- function(numbers, fish,fleet){


  survivors <- vector(mode = 'numeric', length = length(numbers))

  survival  <- exp(-(fish$m + (fleet$eq_f * fleet$sel_at_age)))

  death <-  1 - survival

  survivors[2:fish$max_age] <-
    numbers[1:(fish$max_age - 1)] * survival[1:(fish$max_age - 1)]

  survivors[fish$max_age] <- survivors[fish$max_age] + numbers[fish$max_age] * survival[fish$max_age]

  caught <- (fleet$eq_f * fleet$sel_at_age) / (fish$m + (fleet$eq_f * fleet$sel_at_age)) *  (numbers * death)
  # return(survivors)

  return(list(survivors = survivors, caught = caught))


}
