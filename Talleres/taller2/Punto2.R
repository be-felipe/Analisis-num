library(pracma)
library(Matrix)

A <- matrix(c(-8.1, -7, 6.123, -2, -1, 4, -3, -1, 0, -1, -5, 0.6, -1, 0.33, 6, 1/2), nrow=4, byrow=TRUE)
b <- matrix(c(1.45, 3, 5.12, -4), nrow=4, byrow=TRUE)

diagonal <- function(M){
  M[col(M) != row(M)] <- 0
  return (M)
}

print("Diagonal.")
print(diagonal(A))
print("Triangular Superior.")
print(diagonal(A) - triu(A))
print("Triangular Inferior.")
print(diagonal(A) - tril(A))

#-----------------------Parte b. 

X <- itersolve(A, b, tol=1e-9, method="Gauss-Seidel")
print(X$x)

#-----------------------Parte c. 

x0 <- c(1, 1, 1, 1)
it <- 1

while (it <= 5){
  x1 <- itersolve(A, b, x0, nmax=it, method="Jacobi")
  x2 <- itersolve(A, b, x0, nmax=it + 1, method="Jacobi")
  err <- x2$x - x1$x    
  cat("\nIteracion: ", it, " Error relativo: ", err)
  it <- it + 1
}