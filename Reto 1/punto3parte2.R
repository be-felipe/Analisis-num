
library(PolynomF)

numerosChebyShev = function(intervaloA,intervaloB,listaNumeros)
{
  numerosCalculados= c()
  
  for (i in listaNumeros) 
  {
    numerosCalculados= c(numerosCalculados, 1/2 * (intervaloA + intervaloB)
                         + 1/2 * (intervaloB-intervaloA) 
                         * cos((2*i-1) / (2 * length(listaNumeros))* pi) )
  }
  return (numerosCalculados)
}


n = 3

E = 10e-5

listaCalculada= numerosChebyShev(2^-8,-2^-8,seq(1,n+2))

fxCheby = exp(sin(listaCalculada)-cos(listaCalculada^2))
ns=1
while (ns = ns - 1) 
{
  matrizDeEcuaciones = matrix(nrow = n+2,ncol = n+2)
  
  for (i in 1:(n+2)) 
  {
    matrizDeEcuaciones[i,1] = 1
    
    for (j in 2:(n+1)) 
    {
      matrizDeEcuaciones[i,j] = listaCalculada[i]^(j-1)
    }
    if(i%%2 == 0)
    {
      matrizDeEcuaciones[i,n+2] = -E
    }
    else
    {
      matrizDeEcuaciones[i,n+2] = E
    }
  }
  
  print(matrizDeEcuaciones)
  
  minimax =  solve(matrizDeEcuaciones,fxCheby)
  
  print(minimax)
  
  
  polinomioFinal = polynom(a=minimax[1:(length(minimax)-1)])
  
  print(polinomioFinal)
  
  salto = (2^-8-(-2^-8))/10
  inicio = -pi/64
  errorRelativo = c()
  errorAbsoluto = c()
  for (i in 1:10) {
    
    errorRelativo = c(errorRelativo,abs((exp(sin(listaCalculada)-cos(listaCalculada^2))-polinomioFinal(inicio))/exp(sin(listaCalculada)-cos(listaCalculada^2))))
    errorAbsoluto = c(errorAbsoluto,abs(exp(sin(listaCalculada)-cos(listaCalculada^2))-polinomioFinal(inicio)))
    
    inicio = inicio + salto
    
  }
  
  x = seq(-2^-8,2^-8,length.out = 10)
  
  
  plot(x, exp(sin(x)-cos(x^2)),type = 'b',col="blue", main = "Funcion exp(sin(x)-cos(x^2))" )
  plot(x, polinomioFinal(x),type = 'b',col = "red",main = "Polinomio Generado")
  
  plot(x,errorRelativo,type = 'b',main = "Error Relativo")
  plot(x,errorAbsoluto,type = 'b',main = "Error Absoluto")
  
  
  
}

