
import (matrix)
importh(Math)
library(Matrix)
library(float)
library(PolynomF)

library(pracma)
#Limites de nuestra integra
x1<-0
x2<-0.25
#delta x:
delta<-0.25/4
#Se realizo con 3 rectangulos
#Son los f(x) del problema almacenados en resultado
resultado<-list(0.0625,0.1250,0.2500)
print ( length(resultado))
nuevo<-list()
for (i in 0:length(resultado))
{
 
 temporal<-resultado[i]
# temporal<-sum(temporal,delta)
 nuevo[i]<-temporal
  
}

Resultfinal<-sum(unlist(lapply(nuevo, length)))
Resultfinal<-Resultfinal*delta
print("El resultado utilizando el punto medio de la integral es:")
print(Resultfinal)
