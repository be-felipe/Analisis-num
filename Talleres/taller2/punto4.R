library(pracma)
library(Matrix)

A <- matrix(c( 1, 1/2, 3, 3, 8, 4, 6, 2, 9, 10, 1, 1/4, 1, 12, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14), nrow=5, byrow=TRUE)
b <- matrix(c(4, 6, 0, 4, 8), nrow=5, byrow=TRUE)

multGauss <- function(A, b){
  
  num <- 0
  n <- nrow(A); 
  m <- ncol(A)
  Ab <- cbind(A,b)
  for (k in 1:(n)){
    
    if(k != n){
      fila = which.max( abs(A[k:n,k]) ) + k-1
      
      Ab[c(k, fila), ] = Ab[c(fila, k), ]
      
      if(A[fila,k]==0) stop("La matriz es singular")
      
      for (i in (k+1):n){
        
        Ab[i, ] <- Ab[i, ] - Ab[i, k]/Ab[k,k]*Ab[k, ]
        num <- num + 1
      }  
    }
    for(i in (k):1){
      if(i == k){
        Ab[i, ] <- Ab[i, ]/Ab[i,k]
        num <- num + 1
      }
      else{
        Ab[i, ] <- Ab[i, ] - Ab[i,k]*Ab[k,]
        num <- num + 1
      }
    }
  }
  cat("El numero de multiplicaciones del metodo Gauss es de", num)
}

multGauss(A, b)