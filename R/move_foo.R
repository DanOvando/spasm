#' move function
#'
#' @param numbers
#' @param move_matrix
#'
#' @return
#' @export
#'
move_foo <- function(numbers, move_matrix) {

  moved <- as.numeric(numbers %*% move_matrix)

}