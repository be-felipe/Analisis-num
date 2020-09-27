numeros = c(2,3,5,10,20,30,40,50,60,70)
matriz_triangular <- function (numeros) {
  
  total = c()
  
  for (i in numeros) {
    
    matriz_triangular = matrix (1 , ncol = i, nrow = i)
    aux = 0
    
    for (j in 1:i) {
      for (k in 1:j) {
        
        aux = aux + matriz_triangular[j,k]
      }   
    }
    
    total = c(total , aux)
    aux = 0
  }
  
  return (total)
}

grafica = matriz_triangular (numeros)

plot (numeros , grafica , type = 'l',col = "red",main = "Matriz Triangular ")


# Problema 2.2

numero = c(seq(0,200,by=20))

suma <- function (numero) {
  
  total <- c()
  for (i in numero) {
    
    aux = 0
    
    for (j in 1 : i) {
      
      aux = aux + j^2
    }
    
    total = c (total, aux)
  }
  
  return (total)
}

graficar = suma (numero)
plot (numero , graficar , type = 'l', col ="black",main ="  Sumar numeros cuadrados")



# Problema del Cohete 2.3
f = function()
{
  t = 0
  posicion = function(t){6+(2.13 * (t^2)) - (0.0013 * (t^4))}
  funcion = 6+(2.13 * (t^2)) - (0.0013 * (t^4))
  resultado = 6+2.13*((t+1)^2)-0.0013*((t+1)^4)
  arreglos = c(0)
  
  while(funcion < resultado){
    
    t =t+1
    
    funcion = 6+2.13*(t^2)-0.0013*(t^4)
    
    resultado = 6+2.13*((t+1)^2)-0.0013*((t+1)^4)
    
  }
  plot(posicion, xlim = c(0,10), main = "Trayectoria del Cohete",col ="green")
  
  cat("La altura mayor alcanzada es :",funcion," metros")
  
  
}
f()
