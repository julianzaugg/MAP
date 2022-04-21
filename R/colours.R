

#' Darken a colour
#' @param colour colour name, hexadecimal string or positive integer \code{i} meaning \code{palette()[i]}
#' @param factor value to scale rgb colour vector by division
darken <- function(colour, factor=1.4){
  col <- col2rgb(colour)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#' Lighten a colour
#' @param colour colour name, hexadecimal string or positive integer \code{i} meaning \code{palette()[i]}
#' @param factor value to scale rgb colour vector by multiplication
lighten <- function(colour, factor=1.4){
  col <- col2rgb(colour)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
