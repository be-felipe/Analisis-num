library(pracma)
library(Matrix)

A <- matrix(c(-8.1, -7, 6.123, -2, -1, 4, -3, -1, 0, -1, -5, 0.6, -1, 0.33, 6, 1/2), nrow=4, byrow=TRUE)
b <- matrix(c(1.45, 3, 5.12, -4), nrow=4, byrow=TRUE)

A2 <- matrix(c(4, -1, -1, -1, -1, 4, -1, -1, -1, -1, 4, -1, -1, -1, -1, 4), nrow=4, byrow=TRUE)
b2 <- matrix(c(1, 5, 1.5, -2.33), nrow=4, byrow=TRUE)

diagonal <- function(M){
  M[col(M) != row(M)] <- 0
  return (M)
}

#-----------------------Parte a. Raices del polinomio caracteristico
pol <- charpoly(A, info=TRUE)
raices <- roots(pol$cp)
print("Las raices del polinomio caracteristico de la matriz son:")
print(raices)

#-----------------------Parte c. Matriz de transicion con SOR, Jacobi y Gauss-Seidel

D <- diagonal(A)
U <- D - triu(A)
L <- D - tril(A)
I <- eye(4)
w <- 1

#->SOR
T_sor <- solve((D - w * L)) *((1 - w) * D + w * U)
#El valor optimo de w puede estar en el intervalo [0, 2] ya que dentro de este el metodo
#estara subrelajado y garantiza una mejor convergencia. Si w = 1 el metodo funcionaria 
#igual al metodo de Gauss-Seidel
print("Matriz de transicion con el metodo SOR.")
print(T_sor)

#->Jacobi
T_jac <- -(solve(A)) %*% (L + U)
print("Matriz de transicion con el metodo de Jacobi")
print(T_jac)

#->Gauss-Seidel
T_gs <- (-(solve(A)) %*% U) %*% solve(I + L %*% solve(A))
print("Matriz de transicion con el metodo de Gauss-Seidel")
print(T_gs)

#-----------------------Parte d. Comparacion mejor metodo de convergencia y solucion por defecto

sol_def <- solve(A2, b2)
sol_opt_gs <- itersolve(A2, b2, method="Gauss-Seidel")
sol_opt_jac <- itersolve(A2, b2, method="Jacobi")

print("Solucion del Sistema con solucion por defecto.")
print(sol_def)
cat("Solucion del sistema con solucion G-S con un total de ", sol_opt_gs$iter, " iteraciones.")
print(sol_opt_gs$x)
cat("Solucion del sistema con solucion G-S con un total de ", sol_opt_jac$iter, " iteraciones.")
print(sol_opt_jac$x)