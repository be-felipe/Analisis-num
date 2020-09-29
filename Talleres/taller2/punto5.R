library(pracma)
library(Matrix)

x0 <- c(1, 2, 3)

#--a:

A <- 3
B <- 0

#------------- b:. 

M <- matrix(c(2, 0, -1, B, 2, -1, -1, 1, A), nrow=3, byrow=TRUE)
b <- matrix(c(1, 2, 1), nrow=3, byrow=TRUE)
it <- 1

while(it <= 10){
  x1 <- itersolve(M, b, x0, nmax=it, method="Jacobi")
  cat("\nIteracion: ", it, "\tResultado: ", x1$x)
  it <- it + 1
}