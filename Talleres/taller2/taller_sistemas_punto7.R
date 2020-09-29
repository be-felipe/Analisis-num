
#-----------------------Parte a. Interseccion entre una circunferencia y una recta
interseccion <- function(x){
  n <- length(x)
  f <- rep(NA, n)
  f[1] <- x[1]^2 + x[2]^2 - 1
  f[2] <- x[1] - x[2]
  f
}

x0 <- c(1, 1)
sol = BBsolve(par=x0, fn=interseccion)

plot(sol$par)
plot(interseccion)

#-----------------------Parte b. Analizar el codigo

trigexp = function(x) {
  
  #Tamano del vector que llega por parametro
  n = length(x)
  
  #Se crea un vector F vacio
  F = rep(NA, n)
  
  #Se inserta la primera ecuacion
  F[1] = 3*x[1]^2 + 2*x[2] - 5 + sin(x[1] - x[2]) * sin(x[1] + x[2])
  
  #Se crea una secuencia de 2 a n-1
  tn1 = 2:(n-1)
  
  #Se insertan las siguientes n-1 ecuaciones
  F[tn1] = -x[tn1-1] * exp(x[tn1-1] - x[tn1]) + x[tn1] *
    ( 4 + 3*x[tn1]^2) + 2 * x[tn1 + 1] + sin(x[tn1] -
                                               x[tn1 + 1]) * sin(x[tn1] + x[tn1 + 1]) - 8
  
  #Se inserta la ultima ecuacion
  F[n] = -x[n-1] * exp(x[n-1] - x[n]) + 4*x[n] - 3
  
  #Se retorna F
  F
}

n = 10000
p0 = runif(n) # n initial random starting guesses
sol = BBsolve(par=p0, fn=trigexp) #Se halla la solucion del sistema con un vector n = 1000
sol$par #Se muestra la solucion del sistema