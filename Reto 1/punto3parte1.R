
library(pracma)
library(MASS)

{
  return (exp(sin(x)-cos(x^2)))
}

polinomioCalculado = taylor(f,0,3)


print(as.fractions(polinomioCalculado))

evaluarEnTaylor = function(x,polinomio)
{
  grados = length(polinomio)-1
  total = 0
  for (i in polinomio) 
  {
    total = total + i * (x^grados) 
    grados = grados - 1
  }
  
  return (total)
}


salto = (2^-8-(-2^-8))/10
inicio = -2^-8

for (i in 1:10) {
  cat("Funcion actual ",exp(sin(inicio)-cos(inicio^2)),"\t")
  cat("Valor en Taylor: ",evaluarEnTaylor(inicio,polinomioCalculado),"\n")
  
  inicio = inicio + salto
  
}

