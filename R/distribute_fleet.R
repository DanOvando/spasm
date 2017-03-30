#' \code{distribute_fleet} moves the fishing fleet around
#'
#' @param pop
#' @param effort
#' @param fleet
#'
#' @return a redistributed vector of effort in space
#' @export
#'
#' @examples distribute_fleet(pop, effort, fleet)
distribute_fleet <- function(pop, effort, fleet, num_patches, mpa){


  if (fleet$effort_allocation == 'simple'){
    pop$effort[pop$mpa == F] <- effort / length(unique(pop$patch[pop$mpa == F]))

    # pop$effort <- effort / num_patches

  }

  return(pop$effort)


}