
library(pracma)
library(Matrix)

A <- matrix(c(-8.1, -7, 6.123, -2, -1, 4, -3, -1, 0, -1, -5, 0.6, -1, 0.33, 6, 1/2), nrow=4, byrow=TRUE)

diagonal <- function(M){
  M[col(M) != row(M)] <- 0
  return (M)
}



mat_Ones <- ones(4)
mat_Id <- eye(4)
mat_Zeros <- zeros(4)

print("Matriz de Unos.")
print(mat_Ones)
print("Matriz de Ceros.")
print(mat_Zeros)
print("Matriz Identidad.")
print(mat_Id)

#-Parte b. 
#-> Se comprueba que el determinante de la matriz no sea igual a 0 en primera instansia.
#-> Se comprueba que el determinante de la matriz no sea igual a 0 en primera instansia.

if( det(A) == 0 ){
  print("No hay solucion para el problema. Determinante de A igual a 0. ")
}else{
  I <- mat_Id
  D <- diagonal(A)
  S = upper.tri(A, diag=TRUE)*A
  L = lower.tri(A, diag=TRUE)*A
  Di <- inv(D)
  w <- 1.5 #Factor de Relajacion expresado como w
  
  
  invpri <- inv(I + w*Di*L)
  invpri
  T <- inv(I + w*Di*L) * ((1-w)*(I - w*Di*S))
  print("Matriz de Transicion con el metodo SOR. ")
  print(T)
}