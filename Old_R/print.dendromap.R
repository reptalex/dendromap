#' print dendromap object
#' @export
#' @param x \code{\link{dendromap}} object
#' 
print.dendromap <- function(x){
  print(x$Lineages)
}