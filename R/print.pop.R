#'
#' Generate descriptive summary for objects returned by other functions
#'
#' Generate descriptive summary for objects returned by other functions.
#'
#' @param x the object returned by other functions
#' @param ... ignored arguments
#'
#'
#' @details \code{print()} prints the objects returned by other functions.
#'
#' @return \code{print()} prints the objects returned by other functions.
#'
#'
#'
#'
#' @export

print.pop<-function(x,...){
  print.default(x)
}
