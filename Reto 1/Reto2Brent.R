f <- function(x) x^3-2*x^2+(4*x/3)-(8/27)

pracma::brent(f, 0.5, 1)